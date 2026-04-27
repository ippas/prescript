import logging
import os
import time
import re
import json
import hail as hl
import pandas as pd
from typing import Tuple


class PrescriptionProcessingError(Exception):
    """Raised when prescription filtering fails due to invalid data or logic."""
    pass

class DrugsFiltering:
    """
    Drug prescription filtering utility for UK Biobank General Practice data.

    Version:
        6.1.0

    Overview:
        This class filters and annotates drug prescriptions using predefined code and brand-name dictionaries.
        It supports multiple coding systems and integrates code-based and name-based matching.

    Key features:
        - Dictionary-based filtering using exact drug codes (BNF, Read v2, DM+D).
        - Token-based matching using normalized brand names from drug name prescription field.
        - Paragraph-level BNF codes matching was removed in v6.1.0 in favour of plain drug name matching.
        - Combined output of all matched records, with annotations for identified substances and match mode.

    Expected input Hail table columns:
        `eid`, `source`, `code`, `date`, `system`, `details`

    Each identified prescription is annotated with a list of matched substance names
    and stored in the `substances` field of the resulting table. There are different
    matching modes and the mode in which prescription was matched is in `match_mode`
    field.

    Results are stored in attributes:
        - filtered_prescriptions: combined result of all matched prescriptions
        - unmatched_prescriptions: unmatched rows per system/matching_mode
        - matching_stats: statistics (number) of matched/unmatched records
        - bnf/read2/dmd_code_prescriptions: intermediate results of code matching per system
        - drug_name_prescriptions: remaining prescriptions matched by drug name field
    
    Annotations in output Hail table:
        - `tokenized_drug_name`: Normalized and simplified prescription drug name field.
        - `substances`: List of matched substances for a prescriptions.
        - `match_mode`: Matching mode – `exact_code`, `drug_name`
    """

    PRESCRIPTIONS_COLUMNS = ('idx', 'eid', 'date', 'drug_name', 'quantity', 'tokenized_drug_name', 'substances', 'matched_code', 'match_mode')
    REQUIRED_INPUT_COLUMNS = {'bnf_code': hl.tstr, 'read2_code': hl.tstr, 'dmd_code': hl.tstr, 'drug_name': hl.tstr}

    def __init__(self, raw_input_data: hl.Table, dictionary_dir: str = '../data/codes_lkps', preferred_partitioning: int = 24, logger: logging.Logger = None):
        """
        Initializes a DrugsFiltering instance with input dataset (Hail table), configuration and optional logger.

        Args:
            raw_input_data: Hail Table containing raw prescriptions data.
            dictionary_dir: Directory containing JSON dictionaries (e.g. tokenized_brand_names.json, bnf_codes.json).
            preferred_partitioning (int): Preferred number of partitions for intermediate tables and final table (parallel processing)
            logger: Custom logger. If not provided, a default stderr logger is created.

        Raises:
            TypeError: If any input argument has an incorrect type.
            ValueError: If the input Hail table has an invalid structure.
        """
        if logger is not None and not isinstance(logger, logging.Logger):
            raise TypeError(f'Expected logging.Logger, got {type(logger)} (logger).')
        if not isinstance(raw_input_data, hl.Table):
            raise TypeError(f'Expected hl.Table, got {type(raw_input_data)} (raw_input_data).')
        for column, dtype in self.REQUIRED_INPUT_COLUMNS.items():
            if (column not in raw_input_data.row) or (raw_input_data.row.get(column).dtype != dtype):
                raise ValueError(f'Invalid structure of passed Hail table (raw_input_data.{column}).')
        if not isinstance(dictionary_dir, str):
            raise TypeError(f'Expected str, got {type(dictionary_dir)} (dictionary_dir).')
        if not isinstance(preferred_partitioning, int):
            raise TypeError(f'Expected int, got {type(preferred_partitioning)} (preferred_partitioning).')
        if preferred_partitioning < 1:
            raise ValueError(f'Preferred partition number below 1 (preferred_partitioning).')

        self.logger = logger or self._default_logger(logging.DEBUG)
        
        self.raw_input_data = raw_input_data

        self.dictionary_dir = dictionary_dir

        self.preferred_partitioning = preferred_partitioning

        self.input_data = None
        self.unmatched_prescriptions = None
        self.matching_stats = {}

        self.bnf_code_prescriptions = None
        self.read2_code_prescriptions = None
        self.dmd_code_prescriptions = None
        self.drug_name_prescriptions = None
        self.filtered_prescriptions = None

        self._brand_names_lookup = None
    
    @staticmethod
    def _default_logger(level=logging.INFO) -> logging.Logger:
        """Creates and returns a default stderr logger."""
        logger = logging.getLogger(__name__)
        logger.setLevel(level)
        if not logger.handlers:
            handler = logging.StreamHandler()
            handler.setFormatter(logging.Formatter('%(asctime)s - %(levelname)s - %(message)s'))
            logger.addHandler(handler)
        return logger
    
    def _track_step_execution(method):
        """Decorator that logs duration and handles errors in step execution."""
        def wrapper(*args, **kwargs):
            self = args[0]
            self.logger.info(f'Beginning execution of [{method.__name__}] step...')
            start_time = time.perf_counter()

            try:
                result = method(*args, **kwargs)
            except Exception as e:
                duration = (time.perf_counter() - start_time) / 60.0
                if isinstance(e, PrescriptionProcessingError):
                    self.logger.error(str(e))
                else:
                    msg = " ".join(str(e).splitlines())
                    msg = (msg[:80] + '...') if len(msg) > 80 else msg
                    self.logger.error(f'Unrecognized error: "{msg}" ({type(e).__name__}).')
                self.logger.error(f'Step [{method.__name__}] failed after {duration:.3f} minutes.')
                raise

            duration = (time.perf_counter() - start_time) / 60.0
            self.logger.info(f'Step [{method.__name__}] completed successfully in {duration:.3f} minutes.')

            return result
        
        return wrapper
    
    def _input_data_check(self):
        """Ensures input data has been preprocessed."""
        if self.input_data is None or self.unmatched_prescriptions is None:
            raise PrescriptionProcessingError('Input dataset was not preprocessed (drug name tokenization and code alteration step).')
    
    def _load_dict_json(self, filename: str) -> dict:
        """Loads a JSON dictionary from file and validates its structure."""
        with open(os.path.join(self.dictionary_dir, filename), 'r') as file:
            raw_dict = json.load(file)
        values_count = 0
        if not isinstance(raw_dict, dict):
            raise PrescriptionProcessingError(f'Dictionary loaded from "{filename}" file has invalid structure.')
        for key, values in raw_dict.items():
            if not (isinstance(key, str) and isinstance(values, list) and all(isinstance(item, str) for item in values)):
                raise PrescriptionProcessingError(f'Dictionary loaded from "{filename}" file has invalid structure.')
            values_count += len(values)
        if values_count < 1:
            self.logger.warning(f'Dictionary loaded from file "{filename}" file has no values.')
        self.logger.info(f'Loaded dictionary with {len(raw_dict)} keys and {values_count} values from "{filename}" file ("{self.dictionary_dir}" dir).')
        return raw_dict
    
    def _reverse_dict(self, lookup_dict: dict) -> dict:
        """Reverses a mapping of keys to lists into a mapping of values to keys."""
        reversed_dict = {}
        for key, values in lookup_dict.items():
            for item in values:
                if item not in reversed_dict:
                    reversed_dict[item] = []
                reversed_dict[item].append(key)
        return reversed_dict
    
    def _load_codes_lookup(self, filename: str) -> dict:
        """Loads and reverses code dictionary while checking brand name consistency."""
        codes_dict = self._load_dict_json(filename)
        if any(key not in self.brand_names_lookup for key in codes_dict.keys()):
            raise PrescriptionProcessingError(f'Dictionary loaded from "{filename}" file contains excessive substances (missing in brand names lookup).')
        return self._reverse_dict(codes_dict)
    
    @staticmethod
    def _validate_drug_name_tokenization(tokenized: str) -> bool:
        """Checks if the tokenized drug name matches the expected format."""
        return bool(re.fullmatch(r'^ ((([a-z0-9]*(\d\.\d+)%?[a-z0-9]*)+|[a-z0-9]+) )+$', tokenized))
    
    @property
    def brand_names_lookup(self) -> dict:
        """Returns validated brand name lookup dictionary loaded from file."""
        if self._brand_names_lookup is None:
            brand_names_lkp = self._load_dict_json('tokenized_brand_names.json')
            for brand_names in brand_names_lkp.values():
                for name in brand_names:
                    if not self._validate_drug_name_tokenization(name):
                        raise PrescriptionProcessingError(f'Invalid tokenization of drug brand name "{name}" in dictionary loaded from "brand_names.json" file.')
            self._brand_names_lookup = brand_names_lkp
        return self._brand_names_lookup
    
    @_track_step_execution
    def prepare_input_data(self):
        """
        Filters raw input data to desired coding systems and tokenizes drug names.
        Processed data is repartitioned and cached in memory. Accessible via `input_data` attribute.
        """
        count = self.raw_input_data.count()
        self.matching_stats['raw_input_records'] = count

        self.logger.info(f'{count} records were passed as raw input for prescriptions filtering.')

        raw_data = self.raw_input_data
        
        raw_data = raw_data.repartition(self.preferred_partitioning)

        raw_data = raw_data.annotate(
            tokenized_drug_name = hl.if_else(
                hl.is_defined(raw_data.drug_name),
                hl.rbind(
                    ((' ' + raw_data.drug_name + ' ')
                     .replace(r'(?<=\d),(?=\d{3})', '')
                     .replace(r'(?<!\d)\.(?=\d)', '0.')
                     .replace(r'(?<=\d)\s+%', '%')
                     .replace(r'((?<!\d)\.|\.(?!\d))', ' ')
                     .replace(r'(?<!\d)%', ' ')
                     .replace(r'[\s!"#&\'()*+,\-/:;<=>?@\[\\\]^_`|~]+', ' ')
                     .lower()),
                    lambda tokenized: hl.or_missing(tokenized != ' ', tokenized)
                ),
                hl.missing(hl.tstr)
            )
        )
        
        self.logger.info('Hail drug name tokenization logic for input dataset has been prepared.')

        raw_data = raw_data.annotate(
            bnf_code = hl.if_else(
                hl.is_defined(raw_data.bnf_code) & (hl.len(raw_data.bnf_code) == 14),
                hl.if_else(
                    raw_data.bnf_code.matches(r'^\d{2}\.\d{2}\.\d{2}\.\d{2}\.00$'),
                    raw_data.bnf_code.replace(r'\.', '')[:8],
                    raw_data.bnf_code
                ),
                raw_data.bnf_code)
        )

        raw_data = raw_data.annotate(
            bnf_code = hl.if_else(
                hl.is_defined(raw_data.bnf_code) & (hl.len(raw_data.bnf_code) == 8),
                raw_data.bnf_code[:6],
                raw_data.bnf_code)
        )

        self.logger.info('Hail BNF codes alteration logic for input dataset has been prepared.')

        raw_data = raw_data.annotate(
            read2_code = hl.if_else(
                hl.is_defined(raw_data.read2_code) & (hl.len(raw_data.read2_code) == 7),
                hl.if_else(
                    raw_data.read2_code.matches(r'^[a-z][A-Z0-9a-z.]{3}\.00$'), # slightly faster: read2_data.code[-3:] == '.00'
                    raw_data.read2_code[:5],
                    raw_data.read2_code
                ),
                raw_data.read2_code)
        )

        self.logger.info('Hail Read v2 codes alteration logic for input dataset has been prepared.')

        self.logger.info('Starting drug name tokenization and codes alteration of input dataset...')

        self.input_data = raw_data.cache()
        self.unmatched_prescriptions = self.input_data

        self.logger.info('Input set drug name tokenization and codes alteration completed.')

        count = self.unmatched_prescriptions.count()
        self.logger.info(f'{count} records was prepared for downstream filtering. Results were cached in memory.')
    
    def _hl_dict(self, items: list) -> hl.expr.DictExpression:
        """Converts a list of values into a Hail hash map (used like set)."""
        hashmap = {}
        for c in items:
            hashmap[c] = 1
        return hl.literal(hashmap)
    
    def _build_exact_code_filter(self, system: str, data: hl.Table, codes_lookup: dict) -> Tuple[hl.Table, hl.Table]:
        """Builds filtering and annotating logic to Hail Table using a code-to-substance dictionary."""
        attr_name = f'{system}_code'

        hl_codes_set = self._hl_dict(codes_lookup.keys())

        matched_records = data.filter(hl_codes_set.contains(data[attr_name]))

        hl_codes_lkp = hl.dict(codes_lookup)
        matched_records = matched_records.annotate(
            substances = hl_codes_lkp[matched_records[attr_name]],
            matched_code = hl.struct(system = system, code = matched_records[attr_name]),
            match_mode = f'exact_code'
        )

        unmatched_records = data.filter(~(hl_codes_set.contains(data[attr_name])))

        return matched_records, unmatched_records


    def _filter_exact_code_prescriptions(self, codes_lkp: dict, attr_name: str, system_name: str) -> hl.Table:
        """Filters prescriptions for a given system using exact code lookup."""
        self._input_data_check()

        if getattr(self, f'{attr_name}_code_prescriptions') is not None:
            raise PrescriptionProcessingError(f'Prescriptions have already been filtered by {system_name} codes.')

        self.logger.info(f'Isolated {len(codes_lkp)} {system_name} codes for filtering.')

        count = self.unmatched_prescriptions.count()
        self.logger.debug(f'Remaining {count} records are passed as input for {system_name} codes filtering.')

        code_records, unmatched_records = self._build_exact_code_filter(attr_name, self.unmatched_prescriptions, codes_lkp)

        code_records = code_records.select(*self.PRESCRIPTIONS_COLUMNS)

        self.logger.info(f'Hail filtering and annotating logic for {system_name} system prescriptions has been prepared.')

        self.logger.info(f'Starting execution of {system_name} codes filters...')

        code_records = code_records.cache()
        unmatched_records = unmatched_records.repartition(self.preferred_partitioning).cache()
        system_unmatched_count = unmatched_records.filter(hl.is_defined(unmatched_records[f'{attr_name}_code'])).count()

        setattr(self, f'{attr_name}_code_prescriptions', code_records)
        self.unmatched_prescriptions = unmatched_records

        self.matching_stats[f'{attr_name}_records_not_matched_by_code'] = system_unmatched_count
        self.matching_stats[f'{attr_name}_records_matched_by_code'] = code_records.count()

        self.logger.info(f'{system_name} codes filtering completed. Results were cached in memory.')

        count = unmatched_records.count()
        self.logger.debug(f'Remained {count} records were not matched by {system_name} codes filters. Stored in memory for downstream filtering.')

        count = code_records.count()
        self.logger.info(f'{count} prescriptions were obtained using {system_name} codes filters.')

        return code_records
    
    @_track_step_execution
    def filter_bnf_code_prescriptions(self) -> hl.Table:
        """
        Filters prescriptions using BNF codes based on exact matches.
        Annotates results with matched substances and stores matched and unmatched prescriptions.
        Note: BNF Paragraph-level matching was removed in v6.1.0 in favour of plain drug name matching.
        """
        return self._filter_exact_code_prescriptions(
            self._load_codes_lookup('bnf_codes.json'),
            'bnf',
            'BNF'
        )
    

    @_track_step_execution
    def filter_read2_code_prescriptions(self) -> hl.Table:
        """
        Filters prescriptions using Read v2 codes based on exact matches.
        Annotates results with matched substances and stores matched and unmatched prescriptions.
        """
        return self._filter_exact_code_prescriptions(
            self._load_codes_lookup('read2_codes.json'),
            'read2',
            'Read v2'
        )
    
        
    @_track_step_execution
    def filter_dmd_code_prescriptions(self) -> hl.Table:
        """
        Filters prescriptions using DM+D codes based on exact matches
        Annotates results with matched substances and stores matched and unmatched prescriptions.
        """
        self._filter_exact_code_prescriptions(
            self._load_codes_lookup('dmd_codes.json'),
            'dmd',
            'DM+D'
        )
    
    @_track_step_execution
    def filter_by_drug_name(self) -> hl.Table:
        """
        Filters remaining unmatched prescriptions by tokenized drug names.
        Uses brand name dictionary to annotate substances and stores matched and unmatched prescriptions.
        """
        self._input_data_check()

        if self.drug_name_prescriptions is not None:
            raise PrescriptionProcessingError(f'Prescriptions have already been filtered by drug names.')
              
        count = 0
        for substance_brand_names in self.brand_names_lookup.values():
            count += len(substance_brand_names)
        self.logger.info(f'Isolated {count} drug brand names for name filtering.')

        self.logger.info(f'Preparing unmatched prescriptions for drug name filtering...')
        
        unmatched_code_records = self.unmatched_prescriptions

        name_records = unmatched_code_records.filter(hl.is_defined(unmatched_code_records.tokenized_drug_name))
        name_records = name_records.repartition(self.preferred_partitioning).cache()

        rejected_records = unmatched_code_records.filter(hl.is_missing(unmatched_code_records.tokenized_drug_name))
        
        count = name_records.count()
        self.logger.info(f'Isolated {count} prescriptions for drug name filtering. Data repartitioned and cached in memory.')

        hl_brand_names_lkp = hl.dict(self.brand_names_lookup)
        hl_all_substance = hl.literal(list(self.brand_names_lookup.keys()))

        name_records = name_records.annotate(
            substances = hl_all_substance.filter(
                lambda sub: hl_brand_names_lkp[sub].any(
                    lambda brand_name: name_records.tokenized_drug_name.contains(brand_name))))

        self.logger.info('Hail filtering and annotating logic for drug names prescriptions has been prepared.')
        
        self.logger.info('Starting execution of drug name filters...')

        name_records = name_records.cache()
        
        self.logger.info('Drug name filtering completed. Preparing results and caching in memory...')
        
        matched_name_records = (name_records
                                .filter(hl.len(name_records.substances) > 0)
                                .annotate(
                                    matched_code = hl.missing(hl.tstruct(system=hl.tstr, code=hl.tstr)),
                                    match_mode = 'drug_name')
                                .select(*self.PRESCRIPTIONS_COLUMNS))

        unmatched_name_records = name_records.filter(hl.len(name_records.substances) < 1)
        unmatched_name_records = unmatched_name_records.drop(unmatched_name_records.substances)
        unmatched_records = rejected_records.union(unmatched_name_records)

        matched_name_records = matched_name_records.cache()
        unmatched_records = unmatched_records.cache()

        self.drug_name_prescriptions = matched_name_records
        self.unmatched_prescriptions = unmatched_records

        self.matching_stats['records_not_matched_by_drug_name'] = unmatched_name_records.count()
        self.matching_stats['records_matched_by_drug_name'] = matched_name_records.count()

        self.logger.info('Drug name filtering results are ready and have been cached in memory.')

        count = unmatched_records.count()
        self.logger.debug(f'Remained {count} records were not matched by drug name filters. Stored in memory for downstream filtering.')

        count = matched_name_records.count()
        self.logger.info(f'{count} prescriptions were obtained using drug name filters.')

        return self.drug_name_prescriptions

    @_track_step_execution
    def combine_prescriptions(self) -> hl.Table:
        """
        Joins BNF, Read v2, DM+D f and drug name filtered prescription tables.
        Final result is stored in `filtered_prescriptions` attribute.
        """
        if self.read2_code_prescriptions is None:
            raise PrescriptionProcessingError('Read v2 codes filtering was not executed.')
        if self.bnf_code_prescriptions is None:
            raise PrescriptionProcessingError('BNF codes filtering was not executed.')
        if self.dmd_code_prescriptions is None:
            raise PrescriptionProcessingError('DM+D codes filtering was not executed.')
        if self.drug_name_prescriptions is None:
            raise PrescriptionProcessingError('Drug name filtering was not executed.')

        self.logger.info('Joining all filtered prescriptions...')
        
        joined_prescriptions = (self.read2_code_prescriptions
                                .union(self.bnf_code_prescriptions)
                                .union(self.dmd_code_prescriptions)
                                .union(self.drug_name_prescriptions))
        
        joined_prescriptions = joined_prescriptions.repartition(self.preferred_partitioning).cache()
        
        self.filtered_prescriptions = joined_prescriptions
        self.matching_stats['all_not_matched_records'] = self.unmatched_prescriptions.count()
        self.matching_stats['all_matched_records'] = joined_prescriptions.count()

        self.logger.info('All filtered prescriptions were joined. Results were cached in memory.')

        count = self.filtered_prescriptions.count()
        self.logger.info(f'Total number of filtered prescriptions: {count}.')

        count = self.unmatched_prescriptions.count()
        self.logger.info(f'Total number of remaining not matched prescriptions: {count}.')

        return self.filtered_prescriptions

    
    def run_workflow(self) -> None:
        """
        Runs the full filtering pipeline: prepares input, applies all filters, and merges results.
        The output is saved to the `filtered_prescriptions` attribute.
        """
        self.prepare_input_data()
        self.filter_read2_code_prescriptions()
        self.filter_bnf_code_prescriptions()
        self.filter_dmd_code_prescriptions()
        self.filter_by_drug_name()
        self.combine_prescriptions()




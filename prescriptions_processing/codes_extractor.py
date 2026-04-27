import os
import numpy as np
import pandas as pd
import logging
import requests
import zipfile
import io
import re
import json
from collections import defaultdict
from .dose_annotation import DoseAnnotation

class CodesExtractor:
    def __init__(self, file_prefix=None, output_dir="output", lkps_dir=None):
        """
        Initialize the process with an output directory and optionally a specific LKPs directory, and set up logging.
        
        Args:
            file_prefix (str): Optional prefix for output filenames.
            output_dir (str): Directory to save the output files. Default is "output".
            lkps_dir (str): Directory to save the LKPs files. If None, it will be created in the output directory.
        """
        
        self.logger = logging.getLogger(__name__)
        self.logger.setLevel(logging.INFO)
        handler = logging.StreamHandler()
        handler.setFormatter(logging.Formatter('%(asctime)s - %(levelname)s - %(message)s'))
        self.logger.addHandler(handler)

        self.brand_names_dict = None
        self.omitted_brand_names = None
        self.substances_dict = {}
        self.codes_dicts = {}
        self.doses_and_quantities_dicts = {}

        self.file_prefix = file_prefix
        self.output_dir = output_dir
        os.makedirs(self.output_dir, exist_ok=True)
        self.lkps_dir = lkps_dir if lkps_dir else os.path.join(self.output_dir, "lkps")
        os.makedirs(self.lkps_dir, exist_ok=True)

        self.logger.info(f"Output directory set to: {self.output_dir}")
        self.logger.info(f"LKPs directory set to: {self.lkps_dir}")
        self._ensure_files()

    def _ensure_files(self):
        """
        Check if required files exist, and download them if missing.
        """
        required_files = [
            os.path.join(self.lkps_dir, "bnf_lkp.csv"),
            os.path.join(self.lkps_dir, "dmd_lkp.csv"),
            os.path.join(self.lkps_dir, "read_v2_drugs_lkp.csv")
        ]

        
        if not all(os.path.exists(file) for file in required_files):
            self.logger.info("Required files not found. Downloading...")
            try:
                self._download_files()
            except Exception as e:
                self.logger.error(f"Critical failure when downloading required files: {e}")
                raise
        else:
            self.logger.info("Required files already exist. Skipping download.")


    def _download_files(self):
        """
        Downloads and extracts lookup Excel file from UK Biobank and saves relevant sheets as CSVs.
        """

        url = "https://biobank.ndph.ox.ac.uk/ukb/ukb/auxdata/primarycare_codings.zip"

        response = requests.get(url)
        response.raise_for_status() 

        with zipfile.ZipFile(io.BytesIO(response.content)) as z:
            z.extractall(self.lkps_dir)
            self.logger.info(f"Extracted files to {self.lkps_dir}")
        
        input_excel_file = os.path.join(self.lkps_dir, "all_lkps_maps_v4.xlsx")
        sheets = ["bnf_lkp", "dmd_lkp", "read_v2_drugs_lkp"]
        all_sheets = pd.read_excel(input_excel_file, sheet_name=sheets)
            
        for sheet in sheets:
            csv_file_name = os.path.join(self.lkps_dir, f"{sheet.replace(' ', '_').replace('/', '_')}.csv")
            data = all_sheets[sheet]
            data.to_csv(csv_file_name, index=False)
        
        self.logger.info(f"Saved sheets to files in {self.lkps_dir}")

    def _save_json_dictionary(self, dictionary, filename, dictionary_name, sort=True):
        """
        Sorts a dictionary keys and values and saves it as a JSON file to the output directory.

        Args:
            dictionary (dict): Dictionary to save.
            filename (str): Name of the output JSON file.
            dictionary_name (str): Description for logging.
            sorted (bool): If True, sort values for each key.
        """
        if self.file_prefix:
            filename = f'{self.file_prefix}_{filename}'
        output_file = os.path.join(self.output_dir, filename)
        try:
            if sort == True:
                for key in dictionary.keys():
                    dictionary[key] = sorted(dictionary[key])
            with open(output_file, "w") as json_file:
                json.dump(dictionary, json_file, indent=4, sort_keys=sort)
        except Exception as e:
            self.logger.error(f"Error while saving {dictionary_name} to JSON file: {e}.")
            raise

        self.logger.info(f"Saved {dictionary_name} to {output_file}")
        
    def _save_table(self, df, filename, table_name):
        """
        Saves a pandas DataFrame as a CSV file to the output directory.

        Args:
            df (pd.DataFrame): DataFrame to save.
            filename (str): Name of the output CSV file.
            table_name (str): Description for logging.
        """
        if self.file_prefix:
            filename = f'{self.file_prefix}_{filename}'
        output_file = os.path.join(self.output_dir, filename)
        try:
            df.to_csv(output_file, index=False)
        except Exception as e:
            self.logger.error(f"Error while saving {table_name} to CSV file: {e}.")
            raise

        self.logger.info(f"Saved {table_name} to {output_file}")
    
    @staticmethod
    def _tokenize_drug_name(name: str) -> str:
        """ 
        Preprocesses and tokenizes a drug name string to a simplified lowercase format. 
        """
        name = ' ' + name + ' '
        name = re.sub(r'(?<=\d),(?=\d{3})', '', name)
        name = re.sub(r'(?<!\d)\.(?=\d)', '0.', name)
        name = re.sub(r'(?<=\d)\s+%', '%', name)
        name = re.sub(r'((?<!\d)\.|\.(?!\d))', ' ', name)
        name = re.sub(r'(?<!\d)%', ' ', name)
        name = re.sub(r'[\s!"#&\'()*+,\-/:;<=>?@\[\\\]^_`|~]+', ' ', name)
        if name == ' ':
            return None
        return name.lower()
    
    @staticmethod
    def _compile_refinement_rules(refinement_rules: dict, option: str) -> dict:
        """ 
        Compiles manual refinement rules into regular expressions or tokenized strings. 
        """
        compiled_rules = []
        for substance, rules in refinement_rules.items():
            try:
                regex_text = substance
                substance = re.compile(regex_text)
                if option == 'ignore':
                    for i in range(len(rules)):
                        regex_text = rules[i]
                        rules[i] = re.compile(regex_text)
                elif option == 'add':
                    rules = [CodesExtractor._tokenize_drug_name(rule) for rule in rules]
                    rules = list(filter(lambda rule: rule is not None, rules))
                else:
                    raise ValueError('Invalid refinement option.')
                compiled_rules.append((substance, rules))
            except re.error:
                raise ValueError(f'Manual refinement rules entry "{regex_text}" is not valid regular expression.')
        return compiled_rules
    
    TOKENIZED_INFORMATIVE_SUFFIXES = [
        ' systemic ', ' parent ', ' eye ', ' proprietary preps ', ' gel ', ' inj ', ' soln ', ' buccal ', ' nsl ',
        ' top ', ' oral ', ' mth ', ' rectal ', ' scalp ', ' cap ', ' flushes ', ' crm ', ' sach ', ' syr ', ' tab ',
        ' pharmacia ', ' suppos ', ' oint ', ' spy ', ' susp ', ' liq ', ' loz ', ' inf ', ' pi ', ' pastil ',
        ' anal ', ' ear ', ' vag ', ' blad ', ' inh ', ' dps ']
    
    @classmethod
    def _remove_informative_suffix(cls, brand_name: str) -> str:
        """ 
        Removes popular non-naming, informative suffixes (e.g. tab, cap, inj) from tokenized brand name.
        """
        if not isinstance(brand_name, str):
            raise TypeError(f'Brand name invalid type. Expected str, got {type(brand_name)}.')
        found_suffix = None
        for suffix in cls.TOKENIZED_INFORMATIVE_SUFFIXES:
            if brand_name.endswith(suffix):
                found_suffix = suffix
                break
        if found_suffix is None:
            return brand_name
        brand_name = brand_name[:-(len(found_suffix) - 1)]
        if brand_name.strip() == '':
            raise ValueError('Brand name contains only informative suffixes, no actual name part.')
        return cls._remove_informative_suffix(brand_name)

    def filter_brand_names(self, substances_dict, manual_refinement_rules):
        """
        Extracts, filters, tokenizes and refines brand names from BNF lookup for each substance and its alternative names.
        Any matched but ignored brand names (due to substring matching or manual rules) are recorded in `omitted_brand_names`.

        Args:
            substances_dict (dict): A dictionary where keys are substances and values are lists of their alternative names. The keys must only consist of numbers, letters, and the hyphen (-), with no other special characters.
            manual_refinement_rules (dict): Dict with 'ignore' and/or 'add' rules, containing manual regex-based refinements.
            
        Returns:
             dict: Dictionary mapping each substance to a list of matched and filtered brand names.
        """

        self.brand_names_dict = None
        self.omitted_brand_names = None
        self.substances_dict = substances_dict

        table_path = os.path.join(self.lkps_dir, "bnf_lkp.csv")

        self.logger.info(f"Filtering brand names from {table_path} for given substances.")

        try:
            bnf_table = pd.read_csv(table_path, dtype=str)

            bnf_table = bnf_table.dropna(subset=['BNF_Product', 'BNF_Chemical_Substance'])
            
            bnf_table['BNF_Chemical_Substance'] = bnf_table['BNF_Chemical_Substance'].apply(self._tokenize_drug_name)
            bnf_table = bnf_table.dropna(subset=['BNF_Chemical_Substance']) # data cleaning after tokenization

            brands_names_dict = {}

            for substance, alternatives in substances_dict.items():
                alternatives = set([self._tokenize_drug_name(alt) for alt in ([substance] + alternatives)])
                brands_names_dict[substance] = alternatives

                alternatives = list(alternatives) # alternatives set copying
                
                for alternative in alternatives:

                    matches = bnf_table[
                        bnf_table['BNF_Chemical_Substance'].apply(
                            lambda tokenized_substance: alternative in tokenized_substance
                        )
                    ]

                    if not matches.empty:
                        found_brands = (
                            matches['BNF_Product']
                            .apply(self._tokenize_drug_name)
                            .dropna() # data cleaning after tokenization
                            .apply(self._remove_informative_suffix)
                            .unique()
                            .tolist()
                        )
                        brands_names_dict[substance].update(found_brands)
        
        except Exception as e:
            self.logger.error(f"Error while processing BNF table: {e}")
            raise

        self.logger.info(f"Refining filtered brand names. Removing self-contained (substring) names and applying manual rules.")

        try:
            all_brands = set()
            for brand_names in brands_names_dict.values():
                all_brands.update(brand_names)

            omitted_brand_names = []

            for substance, brand_names in brands_names_dict.items():
                unique_brands = []
                tokenized_substance = self._tokenize_drug_name(substance)
                for brand in brand_names:
                    is_substring = False
                    if brand != tokenized_substance:
                        for other in all_brands:
                            if brand == other:
                                continue
                            elif brand in other:
                                if other not in brand_names:
                                    omitted_brand_names.append((brand, substance, other, 'substring'))
                                    is_substring = True
                                    break
                    if not is_substring:
                        unique_brands.append(brand)
                brands_names_dict[substance] = unique_brands
            
            # Use if not applying substring filtering:
            # for substance in brands_names_dict: brands_names_dict[substance] = list(brands_names_dict[substance]) 

            ignore_rules = self._compile_refinement_rules(manual_refinement_rules.get('ignore', []), 'ignore')
            add_rules = self._compile_refinement_rules(manual_refinement_rules.get('add', []), 'add')
            
            for substance, brand_names in brands_names_dict.items():
                for substance_rule, brand_rules in ignore_rules:
                    if substance_rule.fullmatch(substance) or substance_rule.fullmatch(substance.strip()):
                        filtered_brand_names = []
                        for name in brand_names:
                            if any((rule.fullmatch(name) or rule.fullmatch(name.strip())) for rule in brand_rules):
                                omitted_brand_names.append((name, substance, name, 'manual_refinement'))
                                continue
                            filtered_brand_names.append(name)
                        brand_names = filtered_brand_names
                for substance_rule, extra_names in add_rules:
                    if substance_rule.fullmatch(substance) or substance_rule.fullmatch(substance.strip()):
                        brand_names.extend(extra_names)
                if brand_names is not brands_names_dict[substance]:
                    brands_names_dict[substance] = list(set(brand_names))

                
        except Exception as e:
            self.logger.error(f"Error while filtering obtained brand names: {e}")
            raise

        self._save_json_dictionary(brands_names_dict, 'tokenized_brand_names.json', 'brand names')

        self.brand_names_dict = brands_names_dict
        self.omitted_brand_names = omitted_brand_names

        return brands_names_dict
    
    @staticmethod
    def _clean_code_string(code: str) -> str:
        """ 
        Strips and validates a code string loaded from CSV. 
        """
        if code is None or pd.isna(code):
            return None
        if not isinstance(code, str):
            raise TypeError(f'Expected str, got {type(code)} in code value loaded from CSV lookup.')
        code = code.strip()
        if code == '':
            return None
        return code
    
    @staticmethod
    def _filter_read2_medication_code(code: str) -> bool:
        """ 
        Validates if a Read v2 code is a usable medication (drug) code. 
        """
        if not re.fullmatch(r'^[a-zA-Z0-9.]{5}$', code): # checking proper code form
            raise ValueError(f'Invalid Read v2 code loaded from CSV lookup.')
        if not re.match(r'^[a-z]', code):  # only medications
            return False
        # Read v2 first letter meaning: w - other medical info, y - companies list, p - not-strict medication like inhalers or bandages, x - extended codes
        if code.startswith(('w', 'y')): # filer out other info
            return False
        if code.endswith('...'): # use only substance level codes or above
            return False
        return True

    def extract_codes(self, code_type):
        """
        Extract codes for the given brand names from the specified code type file.

        Args:
            code_type (str): The type of code to extract (e.g., "BNF", "DMD", "READ2").

        Returns:
            dict: A dictionary where keys are substances and values are lists of extracted codes of given type.
        """

        if code_type not in ["BNF", "DMD", "READ2"]:
            raise ValueError(f"Invalid code type: {code_type}. Must be one of ['BNF', 'DMD', 'READ2'].")
        
        code_files = {
            "BNF": ["bnf_lkp.csv", ["code","BNF_Presentation","description","BNF_Chemical_Substance","BNF_Subparagraph","BNF_Paragraph","BNF_Section","BNF_Chapter"]],
            "DMD": ["dmd_lkp.csv", ["code","description"]],
            "READ2": ["read_v2_drugs_lkp.csv", ["code","description", "status"]],
            #"CTV3": ["read_ctv3_lkp.csv", ["code","description", "description_type", "status"]]
        }

        code_file = os.path.join(self.lkps_dir, code_files[code_type][0])
        columns = code_files[code_type][1]

        self.logger.info(f"Extracting codes from {code_file} for given substances.")

        try:
            code_table = pd.read_csv(code_file, dtype=str)
            code_table.columns = columns

            # Input data cleaning:
            code_table['code'] = code_table['code'].apply(self._clean_code_string)
            code_table = code_table.dropna(subset=['code', 'description'])

            code_table['description'] = code_table['description'].apply(self._tokenize_drug_name)
            code_table = code_table.dropna(subset=['description']) # data cleaning after tokenization

            #if code_type == "CTV3" or code_type == 'READ2':
            #    code_table = code_table[code_table['code'].apply(self._filter_read2_ctv3_medication_code)]

            codes_dict = {}

            for substance, brand_names in self.brand_names_dict.items():

                matches = code_table[
                    code_table['description'].apply(
                        lambda tokenized_desc: any((name in tokenized_desc) for name in brand_names)
                    )
                ]

                codes_dict[substance] = matches['code'].unique().tolist()

        except Exception as e:
            self.logger.error(f"Error while processing {code_type} table: {e}")
            raise

        self._save_json_dictionary(codes_dict, f"{code_type.lower()}_codes.json", f'{code_type} codes')

        self.codes_dicts[code_type] = codes_dict

        return codes_dict
    
    
    def _get_code_table(self, code_type):
        code_files = {
            "BNF": ["bnf_lkp.csv", ["code","description","BNF_Product","BNF_Chemical_Substance","BNF_Subparagraph","BNF_Paragraph","BNF_Section","BNF_Chapter"]],
            "DMD": ["dmd_lkp.csv", ["code","description"]],
            "READ2": ["read_v2_drugs_lkp.csv", ["code","description", "status"]],
            #"CTV3": ["read_ctv3_lkp.csv", ["code","description", "description_type", "status"]]
        }
        code_file = os.path.join(self.lkps_dir, code_files[code_type][0])
        columns = code_files[code_type][1]
        code_table = pd.read_csv(code_file, dtype=str)
        code_table.columns = columns
        code_table['code'] = code_table['code'].apply(self._clean_code_string)
        code_table = code_table.dropna(subset=['code', 'description'])
        return code_table
    
    def extract_doses_and_quantities(self, code_type):
        code_table = self._get_code_table(code_type)
        codes_dict = self.codes_dicts[code_type]

        code_to_substances = defaultdict(list)
        for substance, codes in codes_dict.items():
            for code in codes:
                code_to_substances[code].append(substance)

        codes_of_interest = set(code_to_substances.keys())
        filtered_table = code_table[code_table['code'].isin(codes_of_interest)].copy()

        filtered_table['substance'] = filtered_table['code'].map(
            lambda code: sorted(code_to_substances[code])
        )

        filtered_table['norm_description'] = filtered_table['description'].apply(DoseAnnotation.normalize_quantity)
        filtered_table['doses'] = filtered_table['norm_description'].apply(DoseAnnotation.extract_all_doses)
        filtered_table['quantities'] = filtered_table['norm_description'].apply(
            lambda desc: DoseAnnotation.extract_all_quantities(desc, is_second=False)
        )

        filtered_table = filtered_table[['substance', 'code', 'description', 'doses', 'quantities']]
        self.doses_and_quantities_dicts[code_type] = filtered_table
        self._save_table(
            filtered_table[['substance', 'code', 'doses', 'quantities']],
            f"{code_type.lower()}_doses_quantities.csv",
            f"{code_type} doses/quantities table"
        )
        
    def get_all_unique_substance_combinations(self):
        all_combinations = set()
        for df in self.doses_and_quantities_dicts.values():
            if not isinstance(df, pd.DataFrame):
                df = pd.DataFrame(df)
            for subs in df['substance']:
                if len(subs) > 1:
                    combo = tuple(sorted(subs))
                    all_combinations.add(combo)
        return sorted(all_combinations, key=lambda x: (len(x), x))
    
    def get_descriptions_for_given_substance_combinations(self, substance_combinations):
        result = {combo: set() for combo in substance_combinations}
        for df in self.doses_and_quantities_dicts.values():
            if not isinstance(df, pd.DataFrame):
                df = pd.DataFrame(df)
            for _, row in df.iterrows():
                subs = tuple(sorted(row['substance']))
                desc = tuple((row['description'], tuple(row['doses'])))
                if subs in result:
                    result[subs].add(desc)
        return {k: sorted(list(v)) for k, v in result.items()}

    def get_most_common_substance_order(self, substance_combinations):
        combo_doses = defaultdict(lambda: defaultdict(list))
        desc_dict = self.get_descriptions_for_given_substance_combinations(substance_combinations)

        for combo, descriptions in desc_dict.items():
            syn_lists = [self.substances_dict.get(sub, [sub]) for sub in combo]

            for desc_tuple in descriptions:
                description_text, doses = desc_tuple
                if len(doses) != len(combo):
                    continue

                tokens = " " + self._tokenize_drug_name(description_text) + " "
                found = []
                for i, syns in enumerate(syn_lists):
                    for syn in syns:
                        pos = tokens.find(f" {syn} ")
                        if pos != -1:
                            found.append({'pos': pos, 'orig_idx': i})
                            break

                if len(found) != len(combo):
                    continue

                found.sort(key=lambda x: x['pos'])

                for i, found_item in enumerate(found):
                    substance_name = combo[found_item['orig_idx']]
                    try:
                        dose_data = doses[i]
                        dose_val = float(dose_data[0]) if isinstance(dose_data, (list, tuple)) else float(dose_data)
                        if dose_val > 0:
                            combo_doses[combo][substance_name].append(dose_val)
                    except (ValueError, TypeError, IndexError):
                        continue

        result = {}
        for combo, dose_data in combo_doses.items():
            median_doses = []
            for substance_name, dose_list in dose_data.items():
                if dose_list:
                    median = np.median(dose_list)
                    median_doses.append((substance_name, median))

            if len(median_doses) == len(combo):
                median_doses.sort(key=lambda x: x[1])
                ordered_combo = tuple(item[0] for item in median_doses)
                result[combo] = ordered_combo

        result_str_keys = {"|".join(k): list(v) for k, v in result.items()}
        self._save_json_dictionary(result_str_keys, f"substances_order_by_dose.json", f'Order of substances by median dose based on text order', sort=False)

        
    def run_workflow(self, substances_dict, brand_names_refinement_rules):
        """
        Run the entire workflow: download files, filter brand names, and extract codes.
        
        Args:
            substances_dict (dict): A dictionary where keys are substances and values are lists of their alternative names.
            brand_names_refinement_rules (dict): Manual refinement rules for brand name filtering.
        """
        self.filter_brand_names(substances_dict, brand_names_refinement_rules)
        
        code_types = ["BNF", "DMD", "READ2"] # "CTV3"]
        
        for code_type in code_types:
            self.extract_codes(code_type)
        
        for code_type in code_types:
            self.extract_doses_and_quantities(code_type)
    
        substance_combinations = self.get_all_unique_substance_combinations()
        self.get_most_common_substance_order(substance_combinations)
    
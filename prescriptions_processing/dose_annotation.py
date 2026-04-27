import hail as hl
import logging
import time
import re
import os
import pandas as pd
import ast
import json

class DoseAnnotation:    
    def __init__(self, ht: hl.Table, input_dir, logger: logging.Logger = None):
        if logger is not None and not isinstance(logger, logging.Logger):
            raise TypeError(f'Expected logging.Logger, got {type(logger)} (logger).')
        if not isinstance(ht, hl.Table):
            raise TypeError(f'Expected hl.Table, got {type(raw_input_data)} (ht).')
        self.ht = ht
        self.input_dir = input_dir
        self.initial_rows = self.ht.count()
        self.logger = logger or self._default_logger(logging.DEBUG)
        
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
    
    @staticmethod
    def extract_all_doses(text: str) -> list[tuple[str, str]]:
        if text is None:
            return []

        DOSE_UNITS = ["mcg/ml", "mcg", "unit", "%"]
        DOSE_UNIT_REGEX = "|".join([re.escape(u) for u in DOSE_UNITS])
        DOSE_PATTERN = rf'(\d+(\.\d+)?)(\s*)({DOSE_UNIT_REGEX})(?=\s|/|$)'

        return [
            (m.group(1), m.group(4).lower())
            for m in re.finditer(DOSE_PATTERN, text, re.IGNORECASE)
        ]
    
    @staticmethod
    def extract_all_quantities(text: str, is_second: bool) -> list[tuple[str, str]]:
        if text is None:
            return []
        
        if is_second and isinstance(text, str) and re.fullmatch(r'\d+(\.\d+)?', text.strip()):
            return [(text.strip(), "sztuk")]

        results = [
            (m.group(1), "sztuk")
            for m in re.finditer(r'(\d+)\s*sztuk\b', text, re.IGNORECASE)
        ]
        results += [
            (str(int(m.group(1)) * 7), "days")
            for m in re.finditer(r'\b(\d+)\s*weeks?\b', text, re.IGNORECASE)
        ]
        results += [
            (str(int(m.group(1)) * 30), "days")
            for m in re.finditer(r'\b(\d+)\s*months?\b', text, re.IGNORECASE)
        ]

        if not results:
            m = re.search(r'\bx\s*(\d+(\.\d+)?)\b', text, re.IGNORECASE)
            if m:
                return [(m.group(1), "sztuk")]

        return results

    @staticmethod
    def normalize_quantity(text: str) -> str:
        if text is None:
            return ""

        text = text.lower()
        REGEX_DELETE_BRACKETS = r'[\[\]\(\)\{\}\<\>]'
        text = re.sub(REGEX_DELETE_BRACKETS, ' ', text)

        REGEX_MILIGRAMS = r'(?:(?<=\d)|\b|\s)(miligrams|miligram|mg|mgs|milligram|milligrams)\b'
        REGEX_MICROGRAMS = r'(?:(?<=\d)|\b|\s)(micrograms|microgram|Î¼g|Âµg|umg|mcg)\b'
        REGEX_GRAMS = r'(?:(?<=\d)|\b|\s)(gms|g|grams)\b'
        REGEX_MILILITERS = r'(?:(?<=\d)|\b|\s)(milliliter|milliliters|mililiter|mililiters|mls|ml)\b'
        REGEX_UNITS = r'(?:(?<=\d)|\b|\s)(units|unit)\b'
        REGEX_PERCENTAGE = r'(\d+(\.\d+)?)\s*%'
        text = re.sub(REGEX_MILIGRAMS, ' mg ', text, flags=re.IGNORECASE)
        text = re.sub(REGEX_MICROGRAMS, ' mcg ', text, flags=re.IGNORECASE)
        text = re.sub(REGEX_GRAMS, ' g ', text, flags=re.IGNORECASE)
        text = re.sub(REGEX_MILILITERS, ' ml ', text, flags=re.IGNORECASE)
        text = re.sub(REGEX_UNITS, ' unit ', text, flags=re.IGNORECASE)
        text = re.sub(REGEX_PERCENTAGE, r'\1 % ', text, flags=re.IGNORECASE)
    
        units = ["%", "unit", "mg", "mcg", "g", "ml"]
        unit_regex = "|".join(units)

        
        def split_fraction(match):
            left_number = match.group(1)
            right_number = match.group(3)
            unit = match.group(5)
        
            if not unit:
                if right_number == '52':
                    return f"{int(left_number) * 7} days "
                elif right_number == '12':
                    return f"{int(left_number) * 30} days "
                elif right_number == '7':
                    return f"{int(left_number)} days "

            if unit in units:
                return f"{left_number} {unit} / {right_number} {unit} "
            
            return f"{left_number} mg / {right_number} mg "

        FRACTION_WITH_OPTIONAL_UNIT_REGEX = r'(\d+(\.\d+)?)/(\d+(\.\d+)?)(%s)?' % unit_regex
        text = re.sub(
            FRACTION_WITH_OPTIONAL_UNIT_REGEX,
            split_fraction,
            text,
            flags=re.IGNORECASE
        )
        
        def split_mg_ml(match):
            dose_value = float(match.group(1))
            dose_unit = match.group(3).lower()

            return f"{dose_value} {dose_unit}"

        MG_PER_ML_REGEX = r'(\d+(\.\d+)?)\s*(mg|mcg)\s*/\s*(\d+(\.\d+)?)\s*ml'
        text = re.sub(
            MG_PER_ML_REGEX,
            split_mg_ml,
            text,
            flags=re.IGNORECASE
        )

        def convert_to_micrograms(dose_value, unit):
            dose_value = float(dose_value)
            unit = unit.lower()
            if unit == "g":
                return int(dose_value * 1_000_000)
            elif unit == "mg":
                return int(dose_value * 1_000)
            elif unit == "mcg":
                return int(dose_value)
            else:
                return dose_value

        DOSE_TO_MCG_REGEX = r'(\d+(\.\d+)?)(\s*)(g|mg|mcg)\b'
        text = re.sub(
            DOSE_TO_MCG_REGEX,
            lambda m: f"{convert_to_micrograms(m.group(1), m.group(4))} mcg",
            text,
            flags=re.IGNORECASE
        )

        QUANTITIES_UNITS = [
            "suppositories", "suppository", "x piece", "x pieces", "pieces", "piece", "plaster",
            "xtab", "tablest", "tablet", "tablets", "x tablets", "x tablet", "xtablets", "xtablet", "tabs", "tab", "tb", "dispersible tablet", "dispersible tablets", "tabl",
            "x patches", "patches", "patch", "ptch", "doses", "dose", "ampoules", "amp",
            "capsules", "capsule", "cap", "caps", "caplets", "caplet", "ptch", "sach", "lozenges", "sachets", "sachet",            
            "ml", "sztuk"
        ]
        QUANTITIES_UNIT_REGEX = "|".join(QUANTITIES_UNITS)

        TIME_QUANTITES_UNITS = [
            "days", "day", "weeks", "week", "months", "month"
        ]

        def multiply_replacer(m):
            a = int(m.group(1))
            op = m.group(2)
            b = int(m.group(3))
            after = m.group(4) or ""
            if op in ['*', 'x', 'X']:
                result = a * b
        
            after_lower = after.lower()
            if re.match(rf'^\s*({QUANTITIES_UNIT_REGEX})\b', after_lower):
                return f"{result} {after.strip()}"
            else:
                return f"{result} sztuk {after}"
        
        MULTIPLY_QUANTITY_REGEX = r'(\d+)\s*([*xX])\s*(\d+)\b([^\d\w]*[a-zA-Z]+)?'
        text = re.sub(
            MULTIPLY_QUANTITY_REGEX, 
            multiply_replacer, 
            text,
            flags=re.IGNORECASE
        )

        def pack_replacer(m):
            a = int(m.group(1))
            b = int(m.group(2))
            after = (m.group(3) or "").strip()
            after_lower = after.lower()
            if re.match(rf'^\s*({QUANTITIES_UNIT_REGEX})\b', after_lower):
                return f"{a * b} {after}"
            else:
                return f"{a * b} sztuk {after}".strip()
        
        PACK_QUANTITY_REGEX = r'(\d+)\s*pack[s]?\s+of\s*(\d+)\s*([a-zA-Z]+(?:\s*\([^)]+\))?)?'
        text = re.sub(
            PACK_QUANTITY_REGEX,
            pack_replacer,
            text,
            flags=re.IGNORECASE
        )

        QUANTITY_PATTERN = rf'(\b\d+(\.\d+)?)(?:\s*-\s*)?\s*({QUANTITIES_UNIT_REGEX})\b'
        text = re.sub(
            QUANTITY_PATTERN,
            lambda m: f"{m.group(1)} sztuk",
            text,
            flags=re.IGNORECASE
        )

        return text

    def annotate_based_on_description(self):
        unique_pairs = self.ht.aggregate(
            hl.agg.collect_as_set((self.ht.drug_name, self.ht.quantity))
        )

        dose_dict = {}
        quantity_dict = {}
        
        for first, second in unique_pairs:
            norm_first = self.normalize_quantity(first)
            norm_second = self.normalize_quantity(second)
            
            dose_first = self.extract_all_doses(norm_first)
            dose_second = self.extract_all_doses(norm_second)
            dose_merged = list(dose_first)
            set_first = set(dose_first)
            for d in dose_second:
                if d not in set_first:
                    dose_merged.append(d)
            key = (str(first), str(second))        
            dose_dict[key] = dose_merged
           
            quantity_first = self.extract_all_quantities(norm_first, is_second=False)
            quantity_second = self.extract_all_quantities(norm_second, is_second=True)
            quantity = quantity_second
            if quantity_second and len(quantity_second) == 0 and quantity_first and len(quantity_first) > 0:
                quantity = quantity_first
            elif quantity_second and len(quantity_second) > 1:
                quantity = [quantity_second[-1]]
                
            key = (str(first), str(second))
            quantity_dict[key] = quantity

        self.ht = self.ht.annotate(
            doses=hl.literal(dose_dict).get((self.ht.drug_name, self.ht.quantity), hl.empty_array(hl.ttuple(hl.tstr, hl.tstr))),
            quantities=hl.literal(quantity_dict).get((self.ht.drug_name, self.ht.quantity), hl.empty_array(hl.ttuple(hl.tstr, hl.tstr)))
        )
        
        self.ht = self.ht.persist()
        
        doses_annotated_count = self.ht.aggregate(hl.agg.count_where(hl.len(self.ht.doses) > 0))
        quantities_annotated_count = self.ht.aggregate(hl.agg.count_where(hl.len(self.ht.doses) > 0))

        self.logger.info(f"Successfully assigned doses from description to {doses_annotated_count} records ({doses_annotated_count / self.initial_rows:.2%}).")
        self.logger.info(f"Successfully assigned quantities from description to {quantities_annotated_count} records ({quantities_annotated_count / self.initial_rows:.2%}).")
        
        return self.ht
        
        
    def annotate_based_on_code(self):
        code_types = {"BNF": "bnf", "DMD":"dmd", "READ2": "read_2"}

        for code_file, code_hl in code_types.items():
            table_path = os.path.join(self.input_dir, f"{code_file.lower()}_doses_quantities.csv")
            df = pd.read_csv(table_path)
            df['code'] = df['code'].astype(str)
            df['doses'] = df['doses'].apply(ast.literal_eval)
            df['quantities'] = df['quantities'].apply(ast.literal_eval)
            code_to_doses = dict(zip(df['code'], df['doses']))
            code_to_quantities = dict(zip(df['code'], df['quantities']))
            
            self.ht = self.ht.annotate(
                doses=hl.if_else(
                    (self.ht.matched_code.system == code_hl) & ((hl.is_missing(self.ht.doses)) | (hl.len(self.ht.doses) == 0)),
                    hl.literal(code_to_doses).get(hl.str(self.ht.matched_code.code), hl.empty_array(hl.ttuple(hl.tstr, hl.tstr))),
                    self.ht.doses
                ),
                quantities=hl.if_else(
                    (self.ht.matched_code.system == code_hl) & ((hl.is_missing(self.ht.quantities)) | (hl.len(self.ht.quantities) == 0)),
                    hl.literal(code_to_quantities).get(hl.str(self.ht.matched_code.code), hl.empty_array(hl.ttuple(hl.tstr, hl.tstr))),
                    self.ht.quantities
                )
            )
            
        self.ht = self.ht.persist()
        
        final_doses_annotated_count = self.ht.aggregate(hl.agg.count_where(hl.len(self.ht.doses) > 0))
        final_quantities_annotated_count = self.ht.aggregate(hl.agg.count_where(hl.len(self.ht.quantities) > 0))

        self.logger.info(f"Total doses assigned to {final_doses_annotated_count} records ({final_doses_annotated_count / self.initial_rows:.2%}).")
        self.logger.info(f"Total quantities assigned to {final_quantities_annotated_count} records ({final_quantities_annotated_count / self.initial_rows:.2%}).")
        
        return self.ht
        

    def split_quantity_and_dose_columns(self):
        unit_priority = hl.literal(["mcg", "unit", "%"])
        presented_units = hl.set(hl.map(lambda d: d[1], self.ht.doses))
        selected_unit = hl.find(lambda u: presented_units.contains(u), unit_priority)
        filtered = hl.filter(lambda d: d[1] == selected_unit, self.ht.doses)
        dose_values = hl.map(lambda d: hl.int64(hl.float(d[0])), filtered)

        self.ht = self.ht.annotate(
            quantity=hl.if_else(
                hl.is_missing(self.ht.quantities) | (hl.len(self.ht.quantities) == 0),
                hl.struct(is_days=False, value=hl.null(hl.tint64)),
                hl.struct(
                    is_days=self.ht.quantities[0][1] == "days",
                    value=hl.int64(hl.float(self.ht.quantities[0][0]))
                )
            ),
            dose=hl.if_else(
                hl.is_missing(self.ht.doses) | (hl.len(self.ht.doses) == 0),
                hl.struct(unit=hl.null(hl.tstr), values=hl.empty_array(hl.tint64)),
                hl.struct(
                    unit=selected_unit,
                    values=dose_values
                )
            )
        )

        self.ht = self.ht.drop('doses', 'quantities')
        
        return self.ht
    
    
    def load_substance_dict(self, json_dict_path):
        with open(json_dict_path, "r") as f:
            self.combinations_dict = json.load(f)
    
    def annotate_substance_combinations(self):
        self.ht = self.ht.annotate(
            substance_combo=hl.str("|").join(hl.sorted(self.ht.substances)),
            n_substances=hl.len(self.ht.substances)
        )
        
        filtered_ht = self.ht.filter(self.ht.n_substances > 1)

        unique_combos = set(
            filtered_ht.aggregate(hl.agg.collect_as_set(filtered_ht.substance_combo))
        )
        
        dict_keys = set(self.combinations_dict.keys())
        missing = unique_combos - dict_keys
        self.logger.info("Combinations which are in ht but not in the dictionary from UKBB lookups:\n" +
                 "\n\t\t".join(sorted(missing)))
            
    def align_and_explode_substances_and_doses(self):
        self.ht = self.ht.annotate(
            substance_combo=hl.str("|").join(hl.sorted(self.ht.substances))
        )
        
        self.ht = self.ht.annotate(
            aligned_substances=hl.literal(self.combinations_dict).get(self.ht.substance_combo, self.ht.substances)
        )

        self.ht = self.ht.annotate(
            exploded_data = hl.zip(self.ht.aligned_substances, self.ht.dose.values)
                            .map(lambda x: hl.struct(
                                substance=x[0],
                                dose_struct=hl.struct(value=x[1], unit=self.ht.dose.unit)
                            ))
        )


        self.ht = self.ht.explode('exploded_data')

        self.ht = self.ht.annotate(
            substance=self.ht.exploded_data.substance,
            doses=self.ht.exploded_data.dose_struct
        )

        self.ht = self.ht.drop(
            'n_substances',
            'substance_combo',
            'aligned_substances', 
            'exploded_data',     
            'substances',        
            'dose'               
        )
        return self.ht
    
    
    def annotate(self):
        self.logger.info("Starting dose and quantity annotation...")
        self.annotate_based_on_description()
        self.annotate_based_on_code()
        
        self.split_quantity_and_dose_columns()
        substance_dict_path = os.path.join(self.input_dir, 'substances_order_by_dose.json')
        self.load_substance_dict(substance_dict_path)
        self.annotate_substance_combinations()
        self.align_and_explode_substances_and_doses()
        self.logger.info("Dose and quantity annotation finished.")

        return self.ht
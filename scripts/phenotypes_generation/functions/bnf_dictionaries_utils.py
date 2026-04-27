import pandas as pd
from typing import Dict, Any 

def prepare_substance_to_bnf_section_code_dict(filepath: str, substances_column_name = 'base_name', code_column_name = 'bnf_paragraph_code') -> Dict[str, str]:
    '''
    Loads the list of medicines, clears the columns, identifies the BNF section (first 4 digits) 
    and returns a dictionary mapping medicine names to BNF section codes.
    '''
    substances_df = pd.read_csv(filepath, dtype={code_column_name: str, substances_column_name: str} )
    substances_df['base_name'] = substances_df['base_name'].str.split('/')
    substances_df = substances_df.explode('base_name')
    substances_df['base_name'] = substances_df['base_name'].str.strip()
    
    substances_df['bnf_section_code'] = substances_df[code_column_name].str[:4]
    substances_df = substances_df[[substances_column_name, 'bnf_section_code']]
    substances_df = substances_df.drop_duplicates()
    
    grouped_codes = substances_df.groupby('base_name')['bnf_section_code'].apply(lambda x: list(set(x)))
    
    substance_to_code_list_dict = grouped_codes.to_dict()
    
    return substance_to_code_list_dict


def prepare_bnf_section_code_to_short_name_dict(filepath: str, bnf_section_code_column_name = 'bnf_section_code', bnf_short_name_column_name = 'bnf_section_short_name') -> Dict[str, str]:
    
    df = pd.read_csv(filepath, dtype={bnf_section_code_column_name: str, bnf_short_name_column_name: str} )
    grouped_codes = df.groupby(bnf_section_code_column_name)[bnf_short_name_column_name].apply(lambda x: list(set(x)))
    dictionary = grouped_codes.to_dict()
    
    return dictionary
    
    
    
    
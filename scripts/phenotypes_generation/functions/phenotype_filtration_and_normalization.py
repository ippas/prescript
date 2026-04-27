import hail as hl
import pandas as pd
from scipy.stats import boxcox 
from typing import List, Dict, Set, Callable, Tuple

def cap_outliers_mu_plus_n_sigma(
    ht: hl.Table, 
    grouping_col_name: str, 
    phenotype_col_name: str, 
    capped_col_name: str = "capped_phenotype",
    n: float = 8.0
) -> hl.Table:
    """
    Calculates statistics (mean and standard deviation) for the phenotype column 
    within each group (e.g., substance) and trims outliers to the limit 
    mu + n*sigma.
    
    Args:
    ht: Hail table containing the data.
    grouping_col_name: Name of the column by which the data is grouped (e.g. “substance”).
    phenotype_col_name: Name of the column containing the phenotype to be capped.
    capped_col_name: Name of the new column with capped values.

    """

    grouping_col = ht[grouping_col_name]
    phenotype_col = ht[phenotype_col_name]

    substance_stats = ht.group_by(
        grouping_col
    ).aggregate(
        stats=hl.agg.stats(phenotype_col) 
    )
    
    substance_thresholds = substance_stats.annotate(
        mu_plus_n_sigma_limit = substance_stats.stats.mean + (
            n * substance_stats.stats.stdev
        )
    ).persist() 
    
    ht_with_threshold = ht.annotate(
        mu_plus_n_sigma_limit = substance_thresholds[grouping_col].mu_plus_n_sigma_limit
    )
    

    ht_final = ht_with_threshold.annotate(
        **{capped_col_name: hl.min(
            ht_with_threshold[phenotype_col_name],
            ht_with_threshold.mu_plus_n_sigma_limit
        )}
    ).drop('mu_plus_n_sigma_limit') 

    return ht_final.persist()

def delete_outliers_mu_plus_n_sigma(
    ht: hl.Table, 
    grouping_col_name: str, 
    phenotype_col_name: str, 
    n: float = 8.0
) -> hl.Table:
    """
    Calculates statistics (mean and standard deviation) for the phenotype column 
    within each group (e.g., substance) and delete outliers above the limit 
    mu + n*sigma.
    
    Args:
    ht: Hail table containing the data.
    grouping_col_name: Name of the column by which the data is grouped (e.g. “substance”).
    phenotype_col_name: Name of the column containing the phenotype to be capped.

    """

    grouping_col = ht[grouping_col_name]
    phenotype_col = ht[phenotype_col_name]

    substance_stats = ht.group_by(
        grouping_col
    ).aggregate(
        stats=hl.agg.stats(phenotype_col) 
    )
    
    substance_thresholds = substance_stats.annotate(
        mu_plus_n_sigma_limit = substance_stats.stats.mean + (
            n * substance_stats.stats.stdev
        )
    ).persist() 
    
    ht_with_threshold = ht.annotate(
        mu_plus_n_sigma_limit = substance_thresholds[grouping_col].mu_plus_n_sigma_limit
    )
    

    ht_final = ht_with_threshold.filter(
        ht_with_threshold[phenotype_col_name] <= ht_with_threshold.mu_plus_n_sigma_limit
    )
    ht_final = ht_final.drop('mu_plus_n_sigma_limit')

    return ht_final.persist()

def delete_percentile_outliers(
    ht: hl.Table, 
    grouping_col_name: str, 
    phenotype_col_name: str, 
    percentile_val: float = 0.01,  # e.g., 0.01 for the 1st percentile
    n: float = 4.0               
) -> hl.Table:
    """
    Deletes two-sided outliers:
    1. Lower outliers: phenotype < (percentile_val / n)
    2. Upper outliers: phenotype > ((1 - percentile_val) * n)
    Filtering is applied within each group defined by grouping_col_name.

    :param ht: Hail Table (Input Table).
    :param grouping_col_name: Name of the column used for grouping.
    :param phenotype_col_name: Name of the column with values for percentile analysis.
    :param percentile_val: Lower percentile value (e.g., 0.01).
    :param n: Multiplier/divisor used for setting the thresholds.
    :return: Hail Table with two-sided outliers removed.
    """

    # --- Step 1: Count initial records ---
    initial_count = ht.count()

    grouping_col = ht[grouping_col_name]
    phenotype_col = ht[phenotype_col_name]

    # 2. Group and aggregate sorted lists.
    substance_lists = ht.group_by(
        grouping_col
    ).aggregate(
        sorted_phenotype_list=hl.sorted(hl.agg.collect(phenotype_col))
    ).persist()

    # 3. Calculate indices for BOTH lower and upper percentiles (e.g., 1st and 99th).
    upper_percentile_val = 1.0 - percentile_val
    
    substance_indices = substance_lists.annotate(
        list_length = hl.len(substance_lists.sorted_phenotype_list),
        # Lower percentile index (e.g., 1st percentile)
        lower_idx=hl.int(hl.floor(
            (hl.len(substance_lists.sorted_phenotype_list) * percentile_val) 
        )),
        # Upper percentile index (e.g., 99th percentile)
        upper_idx=hl.int(hl.floor(
            (hl.len(substance_lists.sorted_phenotype_list) * upper_percentile_val) 
        ))
    ).persist()

    # 4. Extract both percentile values.
    # Helper for safe index access
    def get_percentile_value(ht, index_col):
        return hl.if_else(
            ht.list_length > 0,
            ht.sorted_phenotype_list[
                hl.min(ht.list_length - 1, 
                       hl.max(0, ht[index_col])
                )
            ],
            hl.missing('float64') 
        )

    substance_percentiles = substance_indices.annotate(
        lower_percentile=get_percentile_value(substance_indices, 'lower_idx'),
        upper_percentile=get_percentile_value(substance_indices, 'upper_idx')
    )

    # 5. Calculate final two-sided limits.
    safe_n = hl.if_else(n == 0.0, hl.float64(hl.missing('float64')), n)

    substance_thresholds = substance_percentiles.annotate(
        # Lower Limit: Percentile / n (Must be GREATER than this to be kept)
        lower_limit=hl.coalesce(substance_percentiles.lower_percentile, hl.float64(0.0)) / safe_n,
        # Upper Limit: Percentile * n (Must be LESS than this to be kept)
        upper_limit=hl.coalesce(substance_percentiles.upper_percentile, hl.float64(0.0)) * n 
    ).select(
        'lower_limit', 'upper_limit' 
    ).persist()
    
    ht_with_threshold = ht.annotate(
        lower_limit = substance_thresholds[grouping_col].lower_limit,
        upper_limit = substance_thresholds[grouping_col].upper_limit
    )
    
    ht_final = ht_with_threshold.filter(
        (ht_with_threshold[phenotype_col_name] >= ht_with_threshold.lower_limit) &
        (ht_with_threshold[phenotype_col_name] <= ht_with_threshold.upper_limit)
    )
    
    ht_final = ht_final.drop('lower_limit', 'upper_limit').persist() 
    
    # --- Step 7: Report results ---
    final_count = ht_final.count()
    removed_count = initial_count - final_count
    
    if initial_count > 0:
        removed_percentage = (removed_count / initial_count) * 100
        print(f"--- Outlier Removal Report ---")
        print(f"Phenotype: {phenotype_col_name}")
        print(f"Parameters: Lower Percentile={percentile_val}, Upper Percentile={upper_percentile_val}, Factor (n)={n}")
        print(f"Total number of records: {initial_count}")
        print(f"Number of records after filtering: {final_count}")
        print(f"Number of removed records: {removed_count}")
        print(f"Percentage of removed records: **{removed_percentage:.4f}%**")
        print(f"---------------------------------------------")
    else:
        print("No records to process.")

    return ht_final


def delete_outliers_above_n_times_percentile(
    ht: hl.Table, 
    grouping_col_name: str, 
    phenotype_col_name: str, 
    percentile_val: float = 0.99, 
    n: float = 4.0               
) -> hl.Table:
    initial_count = ht.count()
    
    grouping_col = ht[grouping_col_name]
    phenotype_col = ht[phenotype_col_name]

    substance_lists = ht.group_by(
        grouping_col
    ).aggregate(
        sorted_phenotype_list=hl.sorted(hl.agg.collect(phenotype_col))
    ).persist()

    substance_thresholds = substance_lists.annotate(
        list_length = hl.len(substance_lists.sorted_phenotype_list),
        percentile_index=hl.int(hl.floor(
            (hl.len(substance_lists.sorted_phenotype_list) * percentile_val) - 1
        ))
    ).persist()

    substance_thresholds = substance_thresholds.annotate(
        percentile_value=hl.if_else(
            substance_thresholds.list_length > 0,
            substance_thresholds.sorted_phenotype_list[hl.max(0, substance_thresholds.percentile_index)],
            hl.missing('float64') 
        )
    )

    substance_thresholds = substance_thresholds.annotate(
        limit=n * hl.coalesce(substance_thresholds.percentile_value, hl.float64(0.0))
    ).select(
        'limit' 
    ).persist()
    
 
    ht_with_threshold = ht.annotate(
        limit = substance_thresholds[grouping_col].limit
    )
    
    ht_final = ht_with_threshold.filter(
        ht_with_threshold[phenotype_col_name] <= ht_with_threshold.limit
    )
    
    ht_final = ht_final.drop('limit')
    
    final_count = ht_final.count()
    
    removed_count = initial_count - final_count
    
    if initial_count > 0:
        removed_percentage = (removed_count / initial_count) * 100
        print(f"--- Outlier Removal Report ---")
        print(f"Total number of records: {initial_count}")
        print(f"Number of records after filtering: {final_count}")
        print(f"Number of removed records: {removed_count}")
        print(f"Percentage of removed records: **{removed_percentage:.4f}%**")
        print(f"---------------------------------------------")
    else:
        print("No records to process.")

    return ht_final.persist()



def pivot_single_phenotype(
    ht: hl.Table, 
    key_col_name: str, 
    pivot_col_name: str, 
    agg_phenotype_col_name: str, 
    col_suffix: str
) -> hl.Table:
    """
    Pivots a long-format table (one row per substance/observation) into a wide-format 
    table (one column per substance/observation type).

    ht: The input Hail Table in long format.
    key_col_name: The name of the column to group by (e.g., 'eid'). This will become the key of the wide-format table.
    pivot_col_name: The name of the key column whose unique values will be expanded into new columns (e.g., 'substance', 'bnf_section_code').
    agg_phenotype_col_name: The name of the phenotype column to aggregate (e.g., 'max_therapy_interval_days').
    col_suffix: Suffix name of the phenotypes columns in output table.
    agg_func: The single aggregation function to use for all phenotypes.
    """
    key_col = ht[key_col_name]
    pivot_col = ht[pivot_col_name]
    phenotype_col = ht[agg_phenotype_col_name]
    
    all_pivot_values = ht.aggregate(hl.agg.collect_as_set(pivot_col))
    
    ht_pivoted = ht.group_by(key_col).aggregate(
        **{f'{sub}_{col_suffix}': hl.agg.filter(
            pivot_col == sub,      
            hl.agg.collect(phenotype_col)
        ) for sub in all_pivot_values}
    )
    
    transmute_exprs = {}
    pivot_columns = [f'{sub}_{col_suffix}' for sub in all_pivot_values]
    
    for col_name in pivot_columns:
        array_col = ht_pivoted[col_name]
        transmute_exprs[col_name] = (
            hl.case()
                .when(hl.len(array_col) > 0, array_col[0]) 
                .default(hl.missing(phenotype_col.dtype))      
        )
    ht_pivoted = ht_pivoted.transmute(**transmute_exprs)
    return ht_pivoted.persist()

def pivot_multiple_phenotypes(
    ht: hl.Table, 
    key_col_name: str, 
    pivot_col_name: str, 
    phenotype_mappings: List[Tuple[str, str]]
) -> hl.Table:
    """
    Pivots a long-format table (one row per substance/observation) into a wide-format 
    table (one column per substance/observation type) for multiple columns.

    ht: The input Hail Table in long format.
    key_col_name: The name of the column to group by (e.g., 'eid'). This will become the key of the wide-format table.
    pivot_col_name: The name of the key column whose unique values will be expanded into new columns (e.g., 'substance', 'bnf_section_code').
    phenotype_mappings: A list of (phenotype_name, suffix) tuples.
    agg_func: The single aggregation function to use for all phenotypes.
    """

    first_phenotype_col, first_suffix = phenotype_mappings[0]
    print(f"Pivoting column '{first_phenotype_col}' with suffix '{first_suffix}'")
    final_ht = pivot_single_phenotype(
        ht, key_col_name, pivot_col_name, first_phenotype_col, first_suffix
    )
    
    final_ht = final_ht.key_by(key_col_name)

    for pheno_col, suffix in phenotype_mappings[1:]:
        print(f"Pivoting column '{pheno_col}' with suffix '{suffix}'")
        current_ht = pivot_single_phenotype(
            ht, key_col_name, pivot_col_name, pheno_col, suffix
        )
        current_ht = current_ht.key_by(key_col_name)
        final_ht = final_ht.join(current_ht, how='outer')
        final_ht = final_ht.persist()
        
    return final_ht.persist()

def phenotype_box_cox_transformation(
    ht: hl.Table, 
    grouping_col_name: str,
    phenotype_col_name: str,
    transformed_col_name: str
) -> hl.Table:
    """
    Applies the Box-Cox power transformation to a phenotype column, 
    calculating the optimal lambda (λ) separately for each defined group.

    Args:
    ht: Hail table containing the data.
    grouping_col_name: Name of the column by which the data is grouped (e.g., 'substance').
    phenotype_col_name: Name of the column containing the phenotype to be transformed.
    transformed_col_name: Name of the transformed column.
    """
    
    working_phenotype_col_name = phenotype_col_name
    
    non_positive_count = ht.aggregate(
        hl.agg.count_where(ht[phenotype_col_name] <= 0)
    )
    
    if non_positive_count > 0:
        C = 1 - ht.aggregate(hl.agg.min(ht[phenotype_col_name])) 
        print(f"WARNING: Non-positive values found. Adding {C:.4f} to the phenotype column.")
        ht = ht.annotate(**{working_phenotype_col_name: ht[phenotype_col_name] + C})
    
    df = ht.key_by().select(
        'eid',
        grouping_col_name,
        working_phenotype_col_name
    ).to_pandas()
    
    def calculate_box_cox_lambda(series):
        data = series.dropna().values.astype(float)
        
        if len(data) < 5: # too few observations
            return None, None
            
        try:
            transformed_data, optimal_lambda = boxcox(data)
            return optimal_lambda, transformed_data
        except ValueError:
            return None, None
        except Exception:
            return None, None

    lambda_map = {}
    
    for name, group in df.groupby(grouping_col_name):
        lambda_value, _ = calculate_box_cox_lambda(group[working_phenotype_col_name])
        lambda_map[name] = lambda_value if lambda_value is not None else 0 

    lambda_hl_dict = hl.literal(lambda_map)

    ht = ht.annotate(
        lambda_val=lambda_hl_dict.get(ht[grouping_col_name])
    )
    
    ht = ht.annotate(**{
        transformed_col_name: hl.if_else(
            ht.lambda_val == 0,
            hl.log(ht[working_phenotype_col_name]),
            ((ht[working_phenotype_col_name] ** ht.lambda_val) - 1) / ht.lambda_val
        )
    })
    
    return ht.drop('lambda_val').persist()
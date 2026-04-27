import numpy as np
import pandas as pd
from scipy import stats

def merge_consecutive_records(records):
    records.sort(key=lambda x: x['date'])
    
    if not records:
        return []

    merged_records = []
    current_group = records[0].copy()
    
    current_group['quantity_value'] = [current_group['quantity_value']]
    current_group['doses_value'] = [current_group['doses_value']]

    for i in range(1, len(records)):
        next_record = records[i].copy()
        date_diff = (next_record['date'] - current_group['date']).days

        if date_diff < 14:
            current_group['interval'] += next_record['interval']
            
            current_group['quantity_value'].append(next_record['quantity_value'])
            current_group['doses_value'].append(next_record['doses_value'])
            
            total_dose_sum = sum(
                q * d for q, d in zip(current_group['quantity_value'], current_group['doses_value'])
            )
            
            current_group['daily_dose'] = total_dose_sum / current_group['interval'] if not pd.isna(current_group['interval']) and current_group['interval'] != 0 else 0
            
        else:
            merged_records.append(current_group)
            current_group = next_record
            current_group['quantity_value'] = [current_group['quantity_value']]
            current_group['doses_value'] = [current_group['doses_value']]

    current_group['sum_quantity'] = sum(current_group['quantity_value'])
    merged_records.append(current_group)
    return merged_records


def detect_dose_changes(daily_dose_series, window_size=8, sd_fraction=2):
    """
    Detects significant changes in a time series of daily doses using a moving median.

    Args:
        daily_dose_series (pd.Series): A time series of daily doses.
        window_size (int, optional): The size of the sliding window for calculating the median. Defaults to 8.
        sd_fraction (int, optional): The standard deviation multiplier to set the threshold for change detection. Defaults to 2.

    Returns:
        np.array: An array of indices where a significant dose change was detected.
    """
    medians = []
    data = daily_dose_series.reset_index(drop=True)

    if data.empty or len(data) < window_size:
        return np.array([])

    for x in range(min((window_size - 1), len(data))):
        medians.append(data.iloc[x])

    for x in range(len(data) - window_size + 1):
        window_median = np.median(data.iloc[x:(x + window_size)])
        medians.append(window_median)

    if len(medians) < 2:
        return np.array([])
    diff_medians = np.diff(medians)

    if len(diff_medians) < window_size:
        return np.array([])
    moving_avg = np.convolve(diff_medians, np.ones(window_size) / window_size, mode='valid')
    sd_diffs = np.std(diff_medians)

    if sd_diffs == 0:
        return np.array([])

    median_threshold = sd_fraction * sd_diffs
    changes = np.abs(diff_medians[window_size - 1:] - moving_avg) > median_threshold
    change_points = np.where(changes)[0] + int(window_size / 2)

    return change_points


def verify_and_filter_changes(df, change_points):
    """
    Verifies detected change points using a t-test to ensure a statistically significant difference.

    Args:
        df (pd.DataFrame): The DataFrame with dose data.
        change_points (np.array): An array of indices from detect_dose_changes.

    Returns:
        np.array: An array of verified change points.
    """
    if len(change_points) == 0:
        return np.array([])

    verified_changes = []
    all_points = np.insert(change_points, 0, -1)
    all_points = np.append(all_points, len(df) - 1)

    for i in range(len(all_points) - 1):
        start_idx = all_points[i] + 1
        end_idx = all_points[i + 1] + 1
        
        try:
            therapy_1 = df.iloc[start_idx:end_idx]
            therapy_2 = df.iloc[end_idx:all_points[i+2]+1]

            if len(therapy_1) > 1 and len(therapy_2) > 1:
                t_stat, p_value = stats.ttest_ind(
                    therapy_1['daily_dose'],
                    therapy_2['daily_dose'],
                    equal_var=False,
                    nan_policy='omit'
                )
                if p_value < 0.05:
                    verified_changes.append(end_idx - 1)
        except (IndexError, ValueError):
            pass

    return np.array(verified_changes)


def process_therapies_with_dose_changes(df, substance, eid, days_to_break=60, sd_f=8, sliding_window=2):
    """
    Analyzes dose data to identify and summarize individual therapies based on treatment breaks and significant dose changes.

    Args:
        df (pd.DataFrame): The input DataFrame containing patient's dose data.
        substance (str): The substance name.
        eid (str): The patient's ID.
        days_to_break (int, optional): The number of days considered a break in treatment. Defaults to 60.
        sd_f (int, optional): The standard deviation fraction for dose change detection. Defaults to 8.
        sliding_window (int, optional): The sliding window size for dose change detection. Defaults to 2.

    Returns:
    pd.DataFrame: A summary DataFrame where each row represents a detected therapy. The DataFrame includes columns
                  with aggregated statistics for each therapy, such as:
                  - 'substance': The drug name.
                  - 'eid': Patient identifier.
                  - 'start_date': The start date of the therapy.
                  - 'end_date': The end date of the therapy.
                  - 'records_in_therapy': The number of records (doses) in the therapy.
                  - 'mean_daily_dose': The average daily dose for the therapy.
                  - 'std_daily_dose': The standard deviation of the daily dose.
                  - 'mse_horizontal_line': Mean squared error of the doses relative to the mean dose.
                  - 'normalized_smae_horizontal_line': Normalized symmetric mean absolute error.
                  - 'normalized_mse_horizontal_line': Normalized mean squared error.
                  - 'therapy_data': A list of dictionaries containing the raw dose data for the therapy.
    """
    if df.empty:
        return pd.DataFrame()

    data = df.copy().reset_index(drop=True)
    if 'sum_quantity' not in data.columns:
        data['sum_quantity'] = 0

    data['is_long_break'] = data['interval'] > (days_to_break + data['sum_quantity'])
    potential_changes = detect_dose_changes(data['daily_dose'], window_size=sliding_window, sd_fraction=sd_f)
    change_indices = verify_and_filter_changes(data, potential_changes)
    
    data['is_dose_change'] = False
    if len(change_indices) > 0:
        data.loc[change_indices, 'is_dose_change'] = True

    data['is_new_therapy'] = data['is_long_break'] | data['is_dose_change']
    data['therapy_id'] = data['is_new_therapy'].cumsum()

    active_therapy_periods = data.copy()
    if active_therapy_periods.empty:
        return pd.DataFrame()

    aggregations = {
        'date': ['min', 'max'],
        'daily_dose': ['mean', 'std', 'count'],
    }
    therapy_summary = active_therapy_periods.groupby('therapy_id').agg(aggregations)
    therapy_summary.columns = ['_'.join(col).strip() for col in therapy_summary.columns.values]
    therapy_summary.rename(columns={
        'date_min': 'start_date',
        'date_max': 'end_date',
        'daily_dose_mean': 'mean_daily_dose',
        'daily_dose_std': 'std_daily_dose',
        'daily_dose_count': 'records_in_therapy'
    }, inplace=True)
    
    therapy_summary['std_daily_dose'].fillna(0, inplace=True)
    
    def get_records(group):
        return group[['daily_dose', 'interval', 'doses_value', 'quantity_value']].to_dict('records')

    therapy_data_struct = active_therapy_periods.groupby('therapy_id').apply(get_records).rename('therapy_data')
    therapy_summary = therapy_summary.join(therapy_data_struct)

    def calculate_normalized_errors(group):
        if group.count() <= 1 or group.mean() == 0:
            return [0.0]

        mean_dose = group.mean()
        normalized_errors = (group - mean_dose) / mean_dose

        return normalized_errors.tolist() 
    
    def calculate_nmse(errors):
        return np.mean(np.square(errors))

    def calculate_nsmae(errors):
        normalized_smae_values = []
        
        for error in errors:
            smae_value = error * np.tanh(error / 2.0) 
            normalized_smae_values.append(smae_value)
        
        return np.mean(normalized_smae_values)


    normalized_errors_data = active_therapy_periods.groupby('therapy_id')['daily_dose'].apply(calculate_normalized_errors).rename('normalized_errors_data')

    therapy_summary = therapy_summary.join(normalized_errors_data)
    therapy_summary['nSMAE'] = therapy_summary['normalized_errors_data'].apply(calculate_nsmae)
    therapy_summary['nMSE'] = therapy_summary['normalized_errors_data'].apply(calculate_nmse)
    
    therapy_summary['substance'] = substance
    therapy_summary['eid'] = eid
    
    if therapy_summary.empty:
        full_period_mean = df['daily_dose'].mean()
        full_period_std = df['daily_dose'].std() if df['daily_dose'].count() > 1 else 0
        normalized_mse = (full_period_std / full_period_mean)**2 if full_period_mean != 0 else np.nan
        errors = (df['daily_dose'] - full_period_mean) / full_period_mean if full_period_mean != 0 else pd.Series(np.nan, index=df.index)
        normalized_smae = (errors * np.tanh(errors / 2)).mean()
        
        single_row_df = pd.DataFrame([{
            'substance': substance,
            'eid': eid,
            'start_date': df['date'].min(),
            'end_date': df['date'].max(),
            'records_in_therapy': df.shape[0],
            'mean_daily_dose': full_period_mean,
            'std_daily_dose': full_period_std,
            'normalized_errors_data': normalized_errors_data,
            'nSMAE': 0.0,
            'nMSE': 0.0,
            'therapy_duration_days': df['sum_quantity'],
            'therapy_data': df.to_dict('records')
        }])
        return single_row_df
    
    therapy_summary['start_date_dt'] = pd.to_datetime(therapy_summary['start_date'], format='%Y-%m-%d') #'2007-11-16'
    therapy_summary['end_date_dt'] = pd.to_datetime(therapy_summary['end_date'], format='%Y-%m-%d')
    therapy_summary['therapy_duration_days'] = (therapy_summary['end_date_dt'] - therapy_summary['start_date_dt']).dt.days
    
    final_cols = [
        'substance', 'eid', 'start_date', 'end_date',
        'records_in_therapy', 'mean_daily_dose', 'std_daily_dose',
        'normalized_errors_data', 'nSMAE', 'nMSE', 'therapy_duration_days', 'therapy_data'
    ]
    
    return therapy_summary.reset_index(drop=True)[final_cols]
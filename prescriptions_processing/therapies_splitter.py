import pandas as pd
import numpy as np
import hail as hl
from scipy import stats
from typing import Optional, Dict, Union

class TherapiesSplitter:
    """
    Groups cleaned prescription records for the same patient and drug into 
    distinct 'therapies' based on time gap and/or statistical dose stability.
    
    Grouping criteria:
    1. Statistically non-significant difference in mean daily dose.
    2. Acceptable time gap between prescriptions X (default 60) _ + quantity.
    """
    
    def __init__(self,
                 ht: hl.Table,
                 time_gap_days: int = 60,
                 use_statistical_dose_check: bool = True,
                 sd_fraction: float = 2.0,
                 window_size: int = 5,
                 batch_size: int = 1000):
        """
        Initializes the splitter with data and parameters.

        Args:
            ht: The input Hail Table with prescriptions.
            time_gap_days: Base time gap (in days) allowed between prescriptions.
            use_statistical_dose_check: If True, uses both time and dose criteria. If False, uses ONLY time criteria.
            sd_fraction: Multiplier for standard deviation used in the statistical dose check (e.g., 2.0 * SD).
            window_size: Number of prescriptions to include in the rolling median/SD calculation for dose check.
        """
        self.ht = ht 
        
        self.time_gap_days = time_gap_days    
        self.use_statistical_dose_check = use_statistical_dose_check
        self.sd_fraction = sd_fraction
        self.window_size = window_size
        self.batch_size = batch_size
        self.all_splits_from_all_batches = []
    
    # --- Helper Methods for Therapy Splitting ---
    
    def _prepare_data(self) -> None:
        self.ht = self.ht.annotate(
            date_struct = hl.struct(
                year = hl.int32(self.ht.date.split('-')[0]),
                month = hl.int32(self.ht.date.split('-')[1]),
                day = hl.int32(self.ht.date.split('-')[2])
            )
        )
        self.ht = self.ht.group_by(self.ht.eid, self.ht.substance).aggregate(
            all_data_collected=hl.agg.collect(
                hl.struct(
                    date_struct=self.ht.date_struct,
                    quantity=self.ht.quantity,
                    doses=self.ht.doses,
                    idx=self.ht.idx
                )
            )
        )

        self.ht = self.ht.annotate(
            all_data_collected_sorted=hl.sorted(
                self.ht.all_data_collected,
                key=lambda s: (s.date_struct.year, s.date_struct.month, s.date_struct.day)
            )
        ).drop('all_data_collected')

        self.ht = self.ht.persist()
        
        self.ht = self.ht.annotate(
            data_with_next=hl.range(0, hl.len(self.ht.all_data_collected_sorted)).map(lambda i:
                hl.struct(
                    date_struct=self.ht.all_data_collected_sorted[i].date_struct,
                    quantity=self.ht.all_data_collected_sorted[i].quantity,
                    doses=self.ht.all_data_collected_sorted[i].doses,
                    idx=self.ht.all_data_collected_sorted[i].idx,
                    next_date_struct=hl.if_else(
                        i < hl.len(self.ht.all_data_collected_sorted) - 1,
                        self.ht.all_data_collected_sorted[i + 1].date_struct,
                        hl.missing(self.ht.all_data_collected_sorted.date_struct.dtype.element_type)
                    )
                )
            )
        ).explode('data_with_next')

        self.ht = self.ht.annotate(
            date_struct=self.ht.data_with_next.date_struct,
            quantity=self.ht.data_with_next.quantity,
            doses=self.ht.data_with_next.doses,
            idx=self.ht.data_with_next.idx,
            interval=hl.if_else(
                hl.is_missing(self.ht.data_with_next.next_date_struct),
                self.ht.data_with_next.quantity.value,
                TherapiesSplitter.days_since_epoch(
                    self.ht.data_with_next.next_date_struct.year,
                    self.ht.data_with_next.next_date_struct.month,
                    self.ht.data_with_next.next_date_struct.day
                ) - TherapiesSplitter.days_since_epoch(
                    self.ht.data_with_next.date_struct.year,
                    self.ht.data_with_next.date_struct.month,
                    self.ht.data_with_next.date_struct.day
                )
            )
        ).drop('data_with_next', 'all_data_collected_sorted')

        self.ht = self.ht.persist()
        
        self.ht = self.ht.group_by(self.ht.eid, self.ht.substance).aggregate(
            all_data_collected=hl.agg.collect(
                hl.struct(
                    date_struct=self.ht.date_struct,
                    quantity=self.ht.quantity,
                    doses=self.ht.doses,
                    idx=self.ht.idx,
                    interval=self.ht.interval
                )
            )
        )

        self.ht = self.ht.annotate(
            all_data_collected_sorted=hl.sorted(
                self.ht.all_data_collected,
                key=lambda s: (s.date_struct.year, s.date_struct.month, s.date_struct.day)
            )
        ).drop('all_data_collected')

        self.ht = self.ht.annotate(
            data_with_prev=hl.range(0, hl.len(self.ht.all_data_collected_sorted)).map(lambda i:
                hl.struct(
                    date_struct=self.ht.all_data_collected_sorted[i].date_struct,
                    quantity=self.ht.all_data_collected_sorted[i].quantity,
                    doses=self.ht.all_data_collected_sorted[i].doses,
                    idx=self.ht.all_data_collected_sorted[i].idx,
                    interval=self.ht.all_data_collected_sorted[i].interval,
                    prev_quantity=hl.if_else(
                        i > 0,
                        self.ht.all_data_collected_sorted[i - 1].quantity,
                        hl.missing(self.ht.all_data_collected_sorted.quantity.dtype.element_type)
                    ),
                    prev_interval=hl.if_else(
                        i > 0,
                        self.ht.all_data_collected_sorted[i - 1].interval,
                        hl.missing(self.ht.all_data_collected_sorted.interval.dtype.element_type)
                    ),
                )
            )
        ).explode('data_with_prev').drop('all_data_collected_sorted')
        
        self.ht = self.ht.persist()
        
        self.ht = self.ht.annotate(
            is_new_therapy = hl.if_else(
                self.ht.data_with_prev.prev_interval > (self.ht.data_with_prev.prev_quantity.value + self.time_gap_days),
                1,
                0
            )
        )
        
        self.ht = self.ht.annotate(
            is_new_therapy = hl.coalesce(self.ht.is_new_therapy, 1),
            date_struct = self.ht.data_with_prev.date_struct,
            quantity = self.ht.data_with_prev.quantity,
            doses = self.ht.data_with_prev.doses,
            idx = self.ht.data_with_prev.idx,
            interval = self.ht.data_with_prev.interval
        ).drop('data_with_prev')

        self.ht = self.ht.persist()
    
        self.ht = self.ht.order_by(
            self.ht.eid, self.ht.substance,
            self.ht.date_struct.year, self.ht.date_struct.month, self.ht.date_struct.day
        )
        
        self.ht = self.ht.annotate(
            tid=hl.scan.count_where(
                hl.bool(self.ht.is_new_therapy)) + self.ht.is_new_therapy
        ).drop('is_new_therapy')
        
        self.ht = self.ht.persist()
        
        return
    
    @staticmethod
    def days_since_epoch(y, m, d):
        days_in_month_array = hl.array([31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31])
        is_leap = (y % 4 == 0) & (y % 100 != 0) | (y % 400 == 0)
        days_in_feb = hl.if_else(is_leap, 29, 28)
        days_from_months = hl.sum(hl.range(0, m - 1).map(lambda i:
            hl.if_else(i == 1, days_in_feb, days_in_month_array[i])
        ))
        days_from_year = hl.sum(hl.range(1970, y).map(lambda year: hl.if_else((year % 4 == 0) & (year % 100 != 0) | (year % 400 == 0), 366, 365)))
        return days_from_year + days_from_months + d
    
    @staticmethod
    def merge_consecutive_records(records):
        if not records:
            return []

        df = pd.DataFrame(records)
        df['date'] = pd.to_datetime(df['date'])
        df.sort_values('date', inplace=True)

        records = df.to_dict('records')

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
                
                if pd.isna(current_group['interval']) or current_group['interval'] == 0:
                    current_group['daily_dose'] = 0
                else:
                    current_group['daily_dose'] = total_dose_sum / current_group['interval']

            else:
                merged_records.append(current_group)
                current_group = next_record
                current_group['quantity_value'] = [current_group['quantity_value']]
                current_group['doses_value'] = [current_group['doses_value']]

        current_group['sum_quantity'] = sum(current_group['quantity_value'])
        merged_records.append(current_group)
        return merged_records
    
    @staticmethod
    def detect_dose_changes(daily_dose_series, window_size, sd_fraction):
        medians = []
        data = daily_dose_series.reset_index(drop=True)

        if data.empty or len(data) < window_size:
            return np.array([])

        for x in range(min((window_size - 1), len(data))):
            medians.append(data.iloc[x])
            
        for x in range(len(data) - window_size + 1):
            window_median = np.nanmedian(data.iloc[x:(x + window_size)])
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
    
    @staticmethod
    def verify_and_filter_changes(df, change_points):
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
    
    @staticmethod
    def find_dose_change_splits(df_patient: pd.DataFrame, tid: int, eid: str, substance: str, sd_f, window_size):    
        records_list = df_patient.to_dict('records')
        merged_data = TherapiesSplitter.merge_consecutive_records(records_list)

        if not merged_data:
            return []

        data = pd.DataFrame(merged_data)
        data['date'] = pd.to_datetime(data['date'])
        data.sort_values('date', inplace=True, ignore_index=True)
        
        data['daily_dose'] = data['daily_dose'].replace(pd.NA, np.nan)

        potential_changes = TherapiesSplitter.detect_dose_changes(data['daily_dose'], window_size=window_size, sd_fraction=sd_f)

        change_indices = TherapiesSplitter.verify_and_filter_changes(data, potential_changes)

        if len(change_indices) == 0:
            return []

        split_dates = data.loc[change_indices, 'date'].tolist()

        results = [(tid, eid, substance, split_date) for split_date in split_dates]

        return results

    def _calculate_daily_doses(self) -> None:
        """
        Calculates the mean daily dose for each prescription and stores it in a new column within self.ht.
        """
        self.ht = self.ht.annotate(
            daily_dose=hl.if_else(
                hl.is_missing(self.ht.quantity) | hl.is_missing(self.ht.doses) | hl.is_missing(self.ht.interval),
                hl.missing(hl.tfloat64),
                hl.float64(self.ht.doses.value) * self.ht.quantity.value / self.ht.interval
            )
        )

        self.ht = self.ht.persist()
        return

    def split(self) -> hl.Table:
        self._prepare_data()
        self._calculate_daily_doses()
        
        ht_grouped_by_therapy = self.ht.group_by('eid', 'substance', 'tid').aggregate(
            n_prescriptions = hl.agg.count(),
            all_data_collected = hl.agg.collect(self.ht.row),
        )
        ht_grouped_by_therapy = ht_grouped_by_therapy.persist()
        
        ht_grouped_by_therapy = ht_grouped_by_therapy.filter(
            ht_grouped_by_therapy.n_prescriptions > self.window_size
        )
        
        num_of_therapies_to_split = ht_grouped_by_therapy.count()
        num_of_batches = (num_of_therapies_to_split + self.batch_size - 1) // self.batch_size
        batch_indices = range(1, num_of_batches + 1)
        
        ht_batched = ht_grouped_by_therapy.annotate(
            row_num = hl.scan.count(), 
        )

        ht_batched = ht_batched.annotate(
            batch_index = (ht_batched.row_num % num_of_batches) + 1
        ).drop('row_num')
        ht_batched = ht_batched.persist()
        
        self.all_splits_from_all_batches = []
        
        for batch_id in batch_indices:
            current_batch_ht = ht_batched.filter(
                ht_batched.batch_index == batch_id
            )
            df_batch = current_batch_ht.explode('all_data_collected').to_pandas()
            df_batch['date'] = pd.to_datetime(
                df_batch['all_data_collected.date_struct.year'].astype(str) + '-' +
                df_batch['all_data_collected.date_struct.month'].astype(str) + '-' +
                df_batch['all_data_collected.date_struct.day'].astype(str)
            )

            df_prepared = df_batch.rename(columns={
                'all_data_collected.quantity.value': 'quantity_value', 
                'all_data_collected.doses.value': 'doses_value',       
                'all_data_collected.interval': 'interval',             
                'all_data_collected.daily_dose': 'daily_dose'          
            })

            df_prepared = df_prepared[[
                'tid',
                'eid', 
                'substance', 
                'date', 
                'quantity_value', 
                'doses_value', 
                'interval', 
                'daily_dose'
            ]].copy()

            df_prepared = df_prepared.dropna(subset=['date'])
            df_prepared.sort_values(by=['tid', 'eid', 'substance', 'date'], inplace=True)

            for (tid_val, eid_val, substance_val), df_group in df_prepared.groupby(['tid', 'eid', 'substance']):
                splits = TherapiesSplitter.find_dose_change_splits(
                    df_group, 
                    tid=tid_val,
                    eid=eid_val, 
                    substance=substance_val, 
                    sd_f=self.sd_fraction, 
                    window_size=self.window_size
                )

                self.all_splits_from_all_batches.extend(splits)
            
        self.all_splits_from_all_batches = pd.DataFrame(self.all_splits_from_all_batches, columns=['tid','eid','substance','date'])
        self.all_splits_from_all_batches['date'] = pd.to_datetime(self.all_splits_from_all_batches['date'])
        self.all_splits_from_all_batches['year'] = self.all_splits_from_all_batches['date'].dt.year
        self.all_splits_from_all_batches['month'] = self.all_splits_from_all_batches['date'].dt.month
        self.all_splits_from_all_batches['day'] = self.all_splits_from_all_batches['date'].dt.day
        self.all_splits_from_all_batches = self.all_splits_from_all_batches.drop('date', axis=1)

        all_splits_from_all_batches_ht = hl.Table.from_pandas(self.all_splits_from_all_batches)
        all_splits_from_all_batches_ht = all_splits_from_all_batches_ht.annotate(
            tid=hl.int64(all_splits_from_all_batches_ht.tid), 
            date_struct=hl.struct(
                year=all_splits_from_all_batches_ht.year,
                month=all_splits_from_all_batches_ht.month,
                day=all_splits_from_all_batches_ht.day
            )
        ).select('tid', 'eid', 'substance', 'date_struct').key_by('tid', 'date_struct')
        all_splits_from_all_batches_ht = all_splits_from_all_batches_ht.persist()
        
        self.ht = self.ht.key_by(
            'tid', 
            date_struct=hl.struct(
                year=self.ht.date_struct.year,
                month=self.ht.date_struct.month,
                day=self.ht.date_struct.day
            )
        )
        self.ht = self.ht.annotate(
            is_new_therapy=hl.is_defined(all_splits_from_all_batches_ht[self.ht.key])
        )
        self.ht = self.ht.order_by(
            self.ht.eid, self.ht.substance,
            self.ht.date_struct.year, self.ht.date_struct.month, self.ht.date_struct.day
        )
        self.ht = self.ht.annotate(
            new_tid_tmp=hl.scan.count_where(
                self.ht.is_new_therapy
            ) + hl.int32(self.ht.is_new_therapy) 
        )
        
        self.ht = self.ht.drop('is_new_therapy')
        self.ht = self.ht.drop('interval')
        
        self.ht = self.ht.rename({
            'new_tid_tmp': 'tid_2',
            'tid': 'tid_1'
        })

        self.ht = self.ht.persist()
        
        key_combination_ht = self.ht.select(
            tid_1=hl.int64(self.ht.tid_1),
            tid_2=hl.int64(self.ht.tid_2)
        ).group_by('tid_1', 'tid_2').aggregate()

        key_combination_ht = key_combination_ht.add_index()
        key_combination_ht = key_combination_ht.annotate(
            final_tid=key_combination_ht.idx + 1
        ).key_by('tid_1', 'tid_2')
        key_combination_ht = key_combination_ht.persist()
        
        self.ht = self.ht.key_by(
            tid_1=hl.int64(self.ht.tid_1),
            tid_2=hl.int64(self.ht.tid_2)
        )

        self.ht = self.ht.annotate(
            tid=key_combination_ht[self.ht.key].final_tid
        )

        self.ht = self.ht.key_by('tid')
        self.ht = self.ht.drop('tid_1', 'tid_2')
        self.ht = self.ht.persist()
        
        return self.ht
            
            

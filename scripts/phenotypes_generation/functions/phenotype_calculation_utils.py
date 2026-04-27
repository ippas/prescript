import hail as hl

def days_since_epoch(y, m, d):
    days_in_month_array = hl.array([31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31])
    is_leap = (y % 4 == 0) & (y % 100 != 0) | (y % 400 == 0)
    days_in_feb = hl.if_else(is_leap, 29, 28)
    days_from_months = hl.sum(hl.range(0, m - 1).map(lambda i:
        hl.if_else(i == 1, days_in_feb, days_in_month_array[i])
    ))
    days_from_year = hl.sum(hl.range(1930, y).map(lambda year: hl.if_else((year % 4 == 0) & (year % 100 != 0) | (year % 400 == 0), 366, 365)))
    return days_from_year + days_from_months + d


def calculate_daily_doses(ht: hl.Table) -> hl.Table:
    ht = ht.annotate(
        total_dose_per_record=ht.quantity.value*ht.doses.value
    )
    ht = ht.group_by('eid', 'substance', 'tid', 'date_struct').aggregate(
        sum_doses = hl.agg.sum(ht.total_dose_per_record),
        sum_quantity = hl.agg.sum(ht.quantity.value),
        prescriptions = hl.agg.collect(
            hl.struct(
                quantity=ht.quantity,                   
                doses=ht.doses,                         
                idx=ht.idx                  
            )
        )
    )
    ht = ht.group_by('eid', 'substance', 'tid').aggregate(
        all_dates_data = hl.agg.collect(
            hl.struct(
                date_struct=ht.date_struct, 
                sum_quantity=ht.sum_quantity,
                sum_doses=ht.sum_doses, 
                prescriptions=ht.prescriptions 
            )
        )
    )
    
    ht = ht.annotate(
        sorted_records = hl.sorted(
            ht.all_dates_data, 
            key=lambda r: (
                r.date_struct.year,
                r.date_struct.month,
                r.date_struct.day
            )
        )
    ).drop('all_dates_data')
    ht = ht.persist()
    
    sample_element = ht.sorted_records[0]
    date_struct_type = sample_element.date_struct.dtype
    missing_date_struct = hl.missing(date_struct_type)
    ht = ht.annotate(
        all_data =hl.range(0, hl.len(ht.sorted_records)).map(lambda i:
             hl.struct(
                prescriptions = ht.sorted_records[i].prescriptions,
                date_struct = ht.sorted_records[i].date_struct,
                sum_quantity = ht.sorted_records[i].sum_quantity,
                sum_doses = ht.sorted_records[i].sum_doses,
                next_date_struct=hl.if_else(
                    i < hl.len(ht.sorted_records) - 1,
                    ht.sorted_records[i + 1].date_struct,
                    missing_date_struct
                )
            )
        )
    ).drop('sorted_records')
    ht = ht.explode('all_data')
    ht = ht.persist()
    
    ht = ht.annotate(
        prescriptions = ht.all_data.prescriptions,
        date_struct = ht.all_data.date_struct,
        sum_quantity = ht.all_data.sum_quantity,
        next_date_struct = ht.all_data.next_date_struct,
        sum_doses = ht.all_data.sum_doses
    ).drop('all_data')
    ht = ht.annotate(
        interval=hl.if_else(
            hl.is_missing(ht.next_date_struct),
            ht.sum_quantity,
            days_since_epoch(
                ht.next_date_struct.year,
                ht.next_date_struct.month,
                ht.next_date_struct.day
            ) - days_since_epoch(
                ht.date_struct.year,
                ht.date_struct.month,
                ht.date_struct.day
            )
        )
    ).drop('next_date_struct')
    ht = ht.annotate(
        daily_dose=ht.sum_doses/ht.interval
    )
    ht = ht.persist()
    
    return ht

def calculate_therapy_doses(ht: hl.Table) -> hl.Table:
    records_to_collect = hl.struct(
        sum_quantity=ht.sum_quantity,
        date_struct=ht.date_struct,
        raw_total_dose=ht.sum_doses,
        interval=ht.interval
    )

    ht_summary = ht.group_by(ht.tid, ht.eid, ht.substance).aggregate(
        all_records = hl.agg.collect(records_to_collect)
    )

    ht_summary = ht_summary.annotate(
        sorted_records = hl.sorted(
            ht_summary.all_records, 
            key=lambda x: hl.tuple([x.date_struct.year, x.date_struct.month, x.date_struct.day])
        ),
    ).drop('all_records')

    ht_summary = ht_summary.annotate(
        min_date_struct = ht_summary.sorted_records[0].date_struct,
        max_date_struct = ht_summary.sorted_records[-1].date_struct,
        duration = hl.sum(ht_summary.sorted_records.interval),
        total_treatment_dose = hl.sum(ht_summary.sorted_records.raw_total_dose)
    )

    ht_summary = ht_summary.annotate(
        daily_dose = hl.if_else(
            ht_summary.duration > 0,
            ht_summary.total_treatment_dose / ht_summary.duration,
            hl.float64(ht_summary.total_treatment_dose) 
        )
    )
    ht_summary = ht_summary.persist()
    
    return ht_summary
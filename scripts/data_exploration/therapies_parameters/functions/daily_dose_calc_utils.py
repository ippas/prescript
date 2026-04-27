import hail as hl

def prepare_and_sort_data(ht: hl.Table) -> hl.Table:
    """
    Groups data by patient and substance, then sorts records chronologically.

    Args:
        ht (hl.Table): The input Hail Table with medication records.

    Returns:
        hl.Table: A grouped Table with a sorted list of all records for each group.
    """
    ht = ht.annotate(
        date_struct = hl.struct(
            year = hl.int32(ht.date.split('-')[0]),
            month = hl.int32(ht.date.split('-')[1]),
            day = hl.int32(ht.date.split('-')[2])
        )
    )
    ht_grouped = ht.group_by(ht.eid, ht.substance).aggregate(
        all_data_collected=hl.agg.collect(
            hl.struct(
                date_struct=ht.date_struct,
                quantity=ht.quantity,
                doses=ht.doses
            )
        )
    )

    ht_grouped = ht_grouped.annotate(
        all_data_collected_sorted=hl.sorted(
            ht_grouped.all_data_collected,
            key=lambda s: (s.date_struct.year, s.date_struct.month, s.date_struct.day)
        )
    ).drop('all_data_collected')
    
    ht_grouped = ht_grouped.persist()
    
    return ht_grouped



def calculate_intervals(ht: hl.Table) -> hl.Table:
    """
    Calculates the time interval in days between consecutive medication records.

    Args:
        ht (hl.Table): A Hail Table from prepare_and_sort_data with a sorted list of records.

    Returns:
        hl.Table: The input Table with an 'interval' field added.
    """
    def days_since_epoch(y, m, d):
        # Implementation from your original script
        days_in_month_array = hl.array([31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31])
        is_leap = (y % 4 == 0) & (y % 100 != 0) | (y % 400 == 0)
        days_in_feb = hl.if_else(is_leap, 29, 28)
        days_from_months = hl.sum(hl.range(0, m - 1).map(lambda i:
            hl.if_else(i == 1, days_in_feb, days_in_month_array[i])
        ))
        days_from_year = hl.sum(hl.range(1970, y).map(lambda year: hl.if_else((year % 4 == 0) & (year % 100 != 0) | (year % 400 == 0), 366, 365)))
        return days_from_year + days_from_months + d

    ht_with_next = ht.annotate(
        data_with_next=hl.range(0, hl.len(ht.all_data_collected_sorted)).map(lambda i:
            hl.struct(
                date_struct=ht.all_data_collected_sorted[i].date_struct,
                quantity=ht.all_data_collected_sorted[i].quantity,
                doses=ht.all_data_collected_sorted[i].doses,
                next_date_struct=hl.if_else(
                    i < hl.len(ht.all_data_collected_sorted) - 1,
                    ht.all_data_collected_sorted[i + 1].date_struct,
                    hl.missing(ht.all_data_collected_sorted.date_struct.dtype.element_type)
                )
            )
        )
    ).explode('data_with_next')

    ht_with_next = ht_with_next.annotate(
        date_struct=ht_with_next.data_with_next.date_struct,
        quantity=ht_with_next.data_with_next.quantity,
        doses=ht_with_next.data_with_next.doses,
        interval=hl.if_else(
            hl.is_missing(ht_with_next.data_with_next.next_date_struct),
            hl.missing(hl.tint32),
            days_since_epoch(
                ht_with_next.data_with_next.next_date_struct.year,
                ht_with_next.data_with_next.next_date_struct.month,
                ht_with_next.data_with_next.next_date_struct.day
            ) - days_since_epoch(
                ht_with_next.data_with_next.date_struct.year,
                ht_with_next.data_with_next.date_struct.month,
                ht_with_next.data_with_next.date_struct.day
            )
        )
    ).drop('data_with_next', 'all_data_collected_sorted')

    ht_with_next = ht_with_next.persist()

    return ht_with_next


def calculate_daily_dose(ht: hl.Table) -> hl.Table:
    """
    Calculates the daily dose for each medication record.

    Args:
        ht (hl.Table): A Hail Table with 'doses' and 'quantity' fields.

    Returns:
        hl.Table: The input Table with a 'daily_dose' field added.
    """
    ht = ht.annotate(
        daily_dose=hl.if_else(
            hl.is_missing(ht.quantity) | hl.is_missing(ht.doses),
            hl.missing(hl.tfloat64),
            hl.float64(ht.doses.value) / ht.quantity.value
        )
    )

    ht = ht.persist()

    return ht
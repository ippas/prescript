import logging
import time
import hail as hl
from typing import Optional, Dict, Union

class DataCleaning:
    """
    Cleans and impute prescription data (hl.Table) based on absurd dose and quantity limits.
    This step should be performed AFTER dose annotation and before TherapiesSplitter (e.g., from DoseAnnotation).
    """

    def __init__(self,
                 ht: hl.Table,
                 daily_dose_limits_in_mg: Optional[Dict[str, Union[int, float]]] = None,
                 impute_values: bool = True,
                 dose_std_dev_threshold: float = 8.0,
                 quantity_std_dev_threshold: float = 8.0,
                 quantity_threshold: int = 1000,
                 logger: logging.Logger = None):
        """
        Initializes the cleaner with data and parameters.

        Args:
            ht: The input Hail Table. Assumes 'dose_per_unit' (or similar) and 'quantity' columns exist.
            daily_dose_limits_in_mg: Optional dict {'drug_name': max_dose_per_tablet} for cleaning.
            dose_std_dev_threshold: Number of standard deviations for outlier cleaning of 'dose' column (if limits not provided).
            quantity_std_dev_threshold: Number of standard deviations for outlier cleaning of 'quantity' column (if limits not provided). 
            quantity_threshold: Maximum acceptable number of pills in the box (absolute limit).
            impute_values: If True, absurd and missing values will be imputed. If False, they and all missing values will be deleted.
        """
        if logger is not None and not isinstance(logger, logging.Logger):
            raise TypeError(f'Expected logging.Logger, got {type(logger)} (logger).')
            
        self.ht = ht
        self.daily_dose_limits = daily_dose_limits_in_mg
        self.dose_std_dev_threshold = dose_std_dev_threshold
        self.quantity_std_dev_threshold = quantity_std_dev_threshold
        self.quantity_threshold = quantity_threshold
        self.impute_values = impute_values
        
        self.initial_count = self.ht.count()
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
        
     # --- Private Methods ---

    def _clean_absurd_doses(self) -> None:
        """
        Cleans the 'dose' column by removing absurdly high values either based on explicit limits or the std_dev_threshold and delete dose values below or equal zero.
        """
        
        self.ht = self.ht.annotate(
            doses=self.ht.doses.annotate(
                value=hl.if_else(
                    self.ht.doses.value <= 0,
                    hl.missing(self.ht.doses.value.dtype),
                    self.ht.doses.value
                ),
                unit=self.ht.doses.unit
            )
        )

        self.ht = self.ht.persist()

        substance_stats = self.ht.group_by(self.ht.substance, self.ht.doses.unit).aggregate(
            dose_stats=hl.agg.stats(self.ht.doses.value),
            count_missing_dose=hl.agg.count_where(~hl.is_missing(self.ht.doses.value))
        )
        substance_stats = substance_stats.persist()

        self.ht = self.ht.annotate(
            stats=substance_stats[(self.ht.substance, self.ht.doses.unit)]
        )
        self.ht = self.ht.persist()

        dose_upper_bound_expr = (
                    self.ht.stats.dose_stats.mean + (self.dose_std_dev_threshold * self.ht.stats.dose_stats.stdev)
        )

        self.ht = self.ht.annotate(
            doses=self.ht.doses.annotate(
                value=hl.if_else(
                    self.ht.doses.value > dose_upper_bound_expr,
                    hl.missing(self.ht.doses.value.dtype),
                    self.ht.doses.value
                )
            ),
        ).drop('stats') 
        self.ht = self.ht.persist()

        
    def _clean_absurd_quantity(self) -> None:
        """
        Cleans the 'quantity' column by removing absurdly high values either based on explicit limit or the std_dev_threshold and delete quantity values below or equal zero.
        """
        self.ht = self.ht.annotate(
            quantity=self.ht.quantity.annotate(
                value=hl.if_else(
                    self.ht.quantity.value <= 0,
                    hl.missing(self.ht.quantity.value.dtype),
                    self.ht.quantity.value
                ),
                is_days=self.ht.quantity.is_days
            )
        )

        self.ht = self.ht.persist()

        substance_stats = self.ht.group_by(self.ht.substance, self.ht.doses.unit).aggregate(
            quantity_stats=hl.agg.stats(self.ht.quantity.value),
            count_missing_quantity=hl.agg.count_where(~hl.is_missing(self.ht.quantity.value))
        )
        substance_stats = substance_stats.persist()

        self.ht = self.ht.annotate(
            stats=substance_stats[(self.ht.substance, self.ht.doses.unit)]
        )
        self.ht = self.ht.persist()

        quantity_upper_bound_expr = (
                    self.ht.stats.quantity_stats.mean + (self.quantity_std_dev_threshold * self.ht.stats.quantity_stats.stdev)
        )

        self.ht = self.ht.annotate(
            quantity=self.ht.quantity.annotate(
                value=hl.if_else(
                    self.ht.quantity.value > self.quantity_std_dev_threshold,
                    hl.missing(self.ht.quantity.value.dtype),
                    self.ht.quantity.value
                )
            ),
        ).drop('stats') 
        self.ht = self.ht.persist()

        self.ht = self.ht.annotate(
            quantity=self.ht.quantity.annotate(
                value=hl.if_else(
                    self.ht.quantity.value > self.quantity_threshold,
                    hl.missing(self.ht.quantity.value.dtype),
                    self.ht.quantity.value
                ),
                is_days=self.ht.quantity.is_days
            )
        )
        self.ht = self.ht.persist()
    
    
    def _impute_missing_values(self) -> None:
        """
        Imputes missing values in 'doses.value' and 'quantity.value' based on the most frequent value for the specific substance/quantity or substance/dose group, but only if the mode represents >=50% of records in that group.
        """
        ht = self.ht
        
        dose_mode_table = ht.group_by(ht.substance, ht.quantity.value, ht.doses.unit).aggregate(
            mode_counts=hl.agg.group_by(ht.doses.value, hl.agg.count()),
            total_records=hl.agg.count()
        )
        
        dose_mode_table = dose_mode_table.annotate(
            most_popular_dose=hl.if_else(
                hl.len(dose_mode_table.mode_counts) > 0,
                hl.sorted(
                    dose_mode_table.mode_counts.items(), 
                    key=lambda x: x[1], 
                    reverse=True
                )[0], 
                hl.missing(hl.ttuple(ht.doses.value.dtype, hl.tint64))
            )
        )
        
        dose_mode_table = dose_mode_table.annotate(
            impute_dose=hl.if_else(
                (dose_mode_table.most_popular_dose[1] / dose_mode_table.total_records) >= 0.5,
                dose_mode_table.most_popular_dose[0],
                hl.null(ht.doses.value.dtype)
            )
        )
        dose_mode_table = dose_mode_table.persist()
        
        quantity_mode_table = ht.group_by(ht.substance, ht.doses.value, ht.doses.unit).aggregate(
            mode_counts=hl.agg.group_by(ht.quantity.value, hl.agg.count()),
            total_records=hl.agg.count()
        )
        
        quantity_mode_table = quantity_mode_table.annotate(
            most_popular_quantity=hl.if_else(
                hl.len(quantity_mode_table.mode_counts) > 0,
                hl.sorted(
                    quantity_mode_table.mode_counts.items(), 
                    key=lambda x: x[1],             
                    reverse=True                        
                )[0],
                hl.missing(hl.ttuple(ht.quantity.value.dtype, hl.tint64))
            )
        )
        quantity_mode_table = quantity_mode_table.annotate(
            impute_quantity=hl.if_else(
                (quantity_mode_table.most_popular_quantity[1] / quantity_mode_table.total_records) >= 0.5,
                quantity_mode_table.most_popular_quantity[0],
                hl.null(ht.quantity.value.dtype)
            )
        )
        quantity_mode_table = quantity_mode_table.persist()
        
        self.ht = self.ht.annotate(
            imputation_lookup=dose_mode_table[self.ht.substance, self.ht.quantity.value, self.ht.doses.unit]
        )
        self.ht = self.ht.annotate(
            doses=self.ht.doses.annotate(
                value=hl.if_else(
                    hl.is_missing(self.ht.doses.value),
                    self.ht.imputation_lookup.impute_dose,
                    self.ht.doses.value
                )
            )
        )
        self.ht = self.ht.drop('imputation_lookup')
        
        self.ht = self.ht.annotate(
            imputation_lookup=quantity_mode_table[self.ht.substance, self.ht.doses.value, self.ht.doses.unit]
        )
        self.ht = self.ht.annotate(
            quantity=self.ht.quantity.annotate(
                value=hl.if_else(
                    hl.is_missing(self.ht.quantity.value),
                    self.ht.imputation_lookup.impute_quantity,
                    self.ht.quantity.value
                )
            )
        )
        self.ht = self.ht.drop('imputation_lookup')
        self.ht = self.ht.persist()

    # --- Public Method ---
        
    def clean_data(self) -> hl.Table:
        """
        Public method to execute all initial data cleaning and optional imputation steps.
        
        Returns:
            The cleaned and potentially imputed Hail Table.
        """
        
        self._clean_absurd_doses()
        self._clean_absurd_quantity()
        
        if self.impute_values:
            self._impute_missing_values()
            
        self.ht = self.ht.filter(
            ~hl.is_missing(self.ht.doses.value) & ~hl.is_missing(self.ht.quantity.value)
        )
        self.ht = self.ht.persist()
        
        self.final_count = self.ht.count()
        deleted_count = self.initial_count - self.final_count

        if self.initial_count > 0:
            deleted_percentage = (deleted_count / self.initial_count) * 100
        else:
            deleted_percentage = 0.0
        
        self.logger.info(
            f"Data cleaning complete. Removed {deleted_count} records "
            f"({deleted_percentage:.2f}%) with missing dose or quantity."
        )

        return self.ht
<div align="center">

  # preSCRIPT
  **Prescription processing framework for pharmacogenomic studies**

</div>

## Overview

preSCRIPT is a framework for large-scale prescription search and annotation in UK Biobank primary-care records. It starts from a user-defined list of target drugs, builds substance-level search dictionaries, finds matching prescriptions in `gp_scripts`, annotates dose and quantity, cleans implausible values, and reconstructs longitudinal therapy episodes.

The framework is designed for pharmacogenomic studies, but the core workflow is drug-list agnostic. You can run it for a completely custom list of substances, for instance:

- your own set of drugs or a single substance
- the full drug panel used in the paper,
- a selected subset of that panel.

The main output of preSCRIPT engine is an annotated Hail table in which prescription rows are assigned to substances, dose and quantity information, and reconstructed into treatment episodes. This table can be used to derive downstream prescription-based phenotypes, such as median daily dose, observed maintenance dose, longest uninterrupted therapy duration, or switches between drugs or classes. Example phenotype derivation workflows are provided in the repository.

## Associated preprint

**preSCRIPT: Large-scale prescription search and annotation engine for pharmacogenomic studies**  
Maria Pieczarka, Paweł Pieńkowski, Paula Konowalska, Sylwia Grubarek, Jacek Hajto, Dzesika Hoinkis, Marcin Piechota, Małgorzata Borczyk and Michał Korostyński.  
Published in **April 2026**.  
DOI: [10.64898/2026.04.28.26351989](https://doi.org/10.64898/2026.04.28.26351989)  
Preprint: [medRxiv preprint](https://www.medrxiv.org/content/10.64898/2026.04.28.26351989v1)

*Journal publication in progress.*

## Requirements

To run the full workflow you need:

- access to UK Biobank primary-care prescription records, especially `gp_scripts` table,
- a Hail/Spark environment,
- the repository files,
- UK Biobank primary-care lookup files,
- a prepared drug list or an already prepared `substances.json`.

The notebooks are written for a DNAnexus / UK Biobank RAP-style environment and write outputs as Hail tables, typically into a DNAnexus-backed database. The same logic can be adapted to another Hail-compatible environment by changing input and output paths.

*This repository does not contain participant-level UK Biobank data, row-level outputs, or restricted lookup files or other restricted data.*

## Key input files

### `data/input/drug_list.csv`

This is the main user-facing drug panel input. It must contain two required columns:

| Column | Meaning |
|---|---|
| `base_name` | Canonical substance name used as the internal substance identifier. |
| `alternative_names` | Comma-separated synonyms, aliases, salt forms, product names, or brand names used to improve search coverage. |

**Combination products** are separated by slashes `/`. Remember to keep the same order of products between columns.

#### Minimal example:

```csv
base_name,alternative_names
amiloride,"Amiloride hydrochloride,Amiloride"
amitriptyline,"Amitriptyline hydrochloride,Amitriptyline"
bendroflumethiazide,"Bendroflumethiazide,bendrofluazide"
codeine/paracetamol,"Codeine phosphate / paracetamol,Co-codamol"
codeine/aspirin,"Codeine phosphate / aspirin,Co-codaprin"
simvastatin/ezetimibe,"Simvastatin and ezetimibe,Simvastatin / Ezetimibe"
venlafaxine,Venlafaxine
warfarin,"Warfarin sodium,Warfarin"
```

#### Important details:

- Values are lowercased and stripped during Step 0.
- The `base_name` itself is automatically added to the list of names for that substance.
- If `base_name` contains `/`, Step 0 treats it as a combination entry and assigns the same alternative names to each component substance. 
- `alternative_names` are separated with commas (not semicolons).
- You can use `alternative_names` strategically to add missing brand names or textual variants that are not well covered by the BNF-derived lookup dictionaries.

Changing `drug_list.csv` is the cleanest way to run preSCRIPT for a different drug panel.

### `data/input/substances.json`

Step 0 converts `drug_list.csv` into `substances.json`. This JSON file is the actual substance dictionary used by Step 1.

Example structure:

```json
{
  "amitriptyline": [
      "amitriptyline"
  ],
  "amlodipine": [
      "amlodipine"
  ],
  "codeine": [
      "co-codamol",
      "co-codaprin",
      "codeine"
  ],
  "warfarin": [
    "warfarin"
  ]
}
```

You do not have to start from `drug_list.csv`. If you already have a valid `substances.json`, you can skip Step 0 and start from Step 1.

### `data/input/brand_names_refinement_rules.json`

This file controls manual refinement of brand-name matching. It is used in Step 1 when preSCRIPT builds brand-name and code dictionaries.

The current structure is:

```json
{
  "ignore": {
    ".+": [
      "numark",
      "boots",
      "lloyds"
    ]
  },
  "add": {
    "citalopram": [
      "cipramil"
    ]
  }
}
```

Use this file as a false-positive control layer:

- `ignore` removes noisy brand-name matches, for instance containing only pharma company name. Regular expressions accepted both as values and keys.
- `add` injects useful missing alternative names manually.
- Keys in `add` or `ignore` matches substance names either directly or via regular expression.
- `".+"` is an example of regular expression, can be used for global ignore rules.

The example rules were tuned for the original drug panel. If you use a new therapeutic area or a heavily modified drug list, review your results in terms of false positive matches or missing prescriptions and alter `brand_names_refinement_rules.json` file accordingly. 

### UK Biobank lookup files

Step 1 uses UK Biobank primary-care lookup files, stored under:

```
data/ukb_lkps/
```

If the required lookup files are missing, the code will download lookups package from UKB and extract the basic primary-care coding lookups. This repository contains lookups for BNF and Read v2 coding systems. However, we do not provide dm+d lookup because of licensing restrictions. But missing dm+d lookup will be automatically downloaded in Step 1.

## Additional documentation

### preSCRIPT technical reference

This README is a practical step-by-step walkthrough for running preSCRIPT. For implementation details, changelog, data-exploration summaries, processing statistics, and additional technical background, see [**`docs/technical_reference.md`**](docs/technical_reference.md).

### UK Biobank Primary Care resources

- Showcase: [https://biobank.ctsu.ox.ac.uk/showcase/label.cgi?id=3000]()
- Primary care data description: [https://biobank.ctsu.ox.ac.uk/showcase/ukb/docs/primary_care_data.pdf]()
- Primary care lookups and detailed documentation: [https://biobank.ctsu.ox.ac.uk/showcase/ukb/auxdata/primarycare_codings.zip]()

## Workflow overview

Run the notebooks in this order:

| Step | Notebook | Main purpose |
|---|---|---|
| Step 0 | `scripts/step_0_prepare_substances_list.ipynb` | Prepare `substances.json` from `drug_list.csv`. |
| Step 1 | `scripts/step_1_extract_brandnames_and_codes.ipynb` | Build brand-name, code, dose, and quantity lookup resources. |
| Step 2 | `scripts/step_2_fetch_dispensed_ukb_data.ipynb` | Materialize the working prescription Hail table from `gp_scripts`. |
| Step 3 | `scripts/step_3_filter_prescriptions.ipynb` | Search prescriptions and retain records matching the target drugs. |
| Step 4 | `scripts/step_4_annotate_dose_and_quantity.ipynb` | Annotate dose and quantity at substance level. |
| Step 5 | `scripts/step_5_clean_and_impute.ipynb` | Clean implausible values and impute missing dose/quantity where possible. |
| Step 6 | `scripts/step_6_split_into_therapies.ipynb` | Split cleaned prescriptions into therapy episodes. |

In the repository we also included resources generated for the drug panel used in the paper, covering cardiovascular and CNS substances. If you only want to test the workflow with this predefined panel, you can skip Steps 0 and 1 and start directly from Step 2, assuming you have access to `gp_scripts`. To run preSCRIPT on a custom drug list, start from Step 0, or from Step 1.

---

## Step 0. Prepare the substance list

**Notebook:** `scripts/step_0_prepare_substances_list.ipynb`

### Input

```text
data/input/drug_list.csv
```

### Output

```text
data/input/substances.json
```

### What happens

This step converts the two-column drug list into the master substance dictionary used by the rest of the framework.

The notebook:

1. reads `drug_list.csv`,
2. keeps `base_name` and `alternative_names`,
3. lowercases and strips values,
4. splits `base_name` on `/`,
5. splits `alternative_names` on commas,
6. assigns combination-product names to each component substance,
7. writes `data/input/substances.json`.

Skip this step if you already prepared `substances.json`.

---

## Step 1. Extracting brand-name and code dictionaries

**Notebook:** `scripts/step_1_extract_brandnames_and_codes.ipynb`  
**Main module:** `prescriptions_processing/codes_extractor.py`

### Input

```text
data/input/substances.json
data/input/brand_names_refinement_rules.json
data/ukb_lkps/
```

### Output

```text
data/codes_lkps/tokenized_brand_names.json
data/codes_lkps/bnf_codes.json
data/codes_lkps/dmd_codes.json
data/codes_lkps/read2_codes.json
data/codes_lkps/bnf_doses_quantities.csv
data/codes_lkps/dmd_doses_quantities.csv
data/codes_lkps/read2_doses_quantities.csv
data/codes_lkps/substances_order_by_dose.json
```

### What happens

This step builds the lookup resources used for prescription search and annotation.

The workflow:

1. loads `substances.json`,
2. loads manual refinement rules,
3. matches substances and aliases against UKB lookup descriptions,
4. builds brand-name dictionary,
5. applies `ignore` and `add` refinement rules,
6. builds BNF, dm+d, and Read v2 code dictionaries,
7. parses lookup descriptions for dose and quantity fallback tables,
8. creates a helper dictionary for assigning doses in combination products.

After this step, preSCRIPT knows which codes and text patterns should be treated as evidence for each target substance.

---

## Step 2. Materialize the source prescription table

**Notebook:** `scripts/step_2_fetch_dispensed_ukb_data.ipynb`

### Input

```text
gp_scripts  # UKB HQL table
```

### Output
```text
prescriptions_db/dispensed_prescriptions.ht  # Hail table
prescriptions_db/dispensed_prescriptions_smp_001.ht  # 1% sample table
prescriptions_db/dispensed_prescriptions_smp_010.ht  # 10% sample table
```

### What happens

This step reads the UK Biobank `gp_scripts` table and writes the local Hail table used by downstream preSCRIPT steps.

The selected fields are:

```text
eid
data_provider -> provider
issue_date    -> date
read_2        -> read2_code
bnf_code
dmd_code
drug_name
quantity
```

The notebook also adds a unique row index:

```text
idx
```

This index makes it easier to trace records across later processing stages.

---

## Step 3. Search prescriptions for the target drugs

**Notebook:** `scripts/step_3_filter_prescriptions.ipynb`  
**Main module:** `prescriptions_processing/drugs_filtering.py`

### Input

```text
prescriptions_db/dispensed_prescriptions.ht  # Hail table
data/codes_lkps/   # Step 1 dictionaries
```

### Output

```text
prescriptions_db/filtered_prescriptions_v6.2.0.ht  # Hail table
```

### What happens

This is the main prescription search step. It identifies rows from the source prescription table that match your target drug panel.

The workflow:

1. reads `dispensed_prescriptions.ht`,
2. normalizes drug names and codes,
3. matches prescriptions by Read v2 code,
4. matches prescriptions by BNF code,
5. matches prescriptions by dm+d code,
6. searches still-unmatched rows using tokenized `drug_name`,
7. writes the matched prescription table.

Main added fields include:

```text
tokenized_drug_name
substances
matched_code
match_mode
```

The matching strategy intentionally uses code-based matching first and free-text rescue afterward. This gives specificity from exact codes while still recovering prescriptions where useful information is present only in `drug_name`.

This is where the contents of `drug_list.csv` and `substances.json` strongly affect recall. If a brand name or synonym is missing from the input dictionary and not present in UKB lookups, preSCRIPT may not find those prescriptions.

---

## Step 4. Annotate dose and quantity

**Notebook:** `scripts/step_4_annotate_dose_and_quantity.ipynb`  
**Main module:** `prescriptions_processing/dose_annotation.py`

### Input

```text
prescriptions_db/filtered_prescriptions_v6.2.0.ht  # Hail table
data/codes_lkps/  # Step 1 dictionaries
```

### Output

```text
prescriptions_db/filtered_prescriptions_with_doses_v6.2.0.ht  # Hail table
```

### What happens

This step adds structured dose and quantity annotations.

The workflow:

1. parses dose and quantity from prescription text fields,
2. normalizes units,
3. extracts candidate dose and quantity values,
4. uses code-based fallback lookup tables when text parsing is insufficient,
5. aligns dose values to substances in combination products,
6. converts matched prescriptions into substance-level rows where needed.

After this step, matched prescriptions are closer to the structure needed for longitudinal exposure analysis: one row should correspond to one substance-level prescription event with interpretable dose and quantity fields where possible.

---

## Step 5. Clean and impute dose/quantity values

**Notebook:** `scripts/step_5_clean_and_impute.ipynb`  
**Main module:** `prescriptions_processing/data_cleaning.py`

### Input

```text
prescriptions_db/filtered_prescriptions_with_doses_v6.2.0.ht  # Hail table
```

### Output

```text
prescriptions_db/cleaned_prescriptions_with_doses_v6.2.0.ht  # Hail table
```

### What happens

This step cleans the dose-annotated prescription table before therapy reconstruction.

The notebook initializes `DataCleaning` with the settings listed below. These defaults were selected empirically based on exploratory statistical analyses. We encourage users to adjust them to their own drug panels, datasets, and research questions.

```python
dose_std_dev_threshold=8
quantity_std_dev_threshold=10000
quantity_threshold=1000
impute_values=True
```

The cleaning logic:

1. removes non-positive dose values,
2. removes non-positive quantity values,
3. removes or blanks extreme outliers,
4. imputes missing dose or quantity values where the within-substance evidence is strong enough,
5. drops records that remain unresolved after cleaning and imputation,
6. writes the cleaned Hail table.

The output is the cleaned prescription-level table used directly by Step 6.

---

## Step 6. Split prescriptions into therapies

**Notebook:** `scripts/step_6_split_into_therapies.ipynb`  
**Main module:** `prescriptions_processing/therapies_splitter.py`

### Input

```text
prescriptions_db/cleaned_prescriptions_with_doses_v6.2.0.ht  # Hail table
```

### Output

```text
prescriptions_db/cleaned_prescriptions_splited_to_therapies_v6.2.0.ht  # Hail table
```

### What happens

This step reconstructs therapy episodes for each participant and substance.

The notebook initializes `TherapiesSplitter` with the parameters listed below. These defaults were selected empirically based on exploratory statistical analyses. We encourage users to adjust them to their own drug panels, datasets, and research questions.

```python
time_gap_days=60
use_statistical_dose_check=True
sd_fraction=2.0
window_size=5
batch_size=10000
```

The splitter:

1. groups prescriptions by participant and substance,
2. orders prescriptions by date,
3. estimates inter-prescription intervals,
4. starts a new therapy when the medication gap is too large,
5. calculates daily dose where possible,
6. detects dose shifts using a sliding window,
7. assigns final therapy identifiers.

The final table contains cleaned, dose-annotated, substance-level prescriptions assigned to reconstructed therapy episodes.

---

## Main outputs

After running all steps, the key Hail tables are:

```text
prescriptions_db/dispensed_prescriptions.ht
prescriptions_db/filtered_prescriptions_v6.2.0.ht
prescriptions_db/filtered_prescriptions_with_doses_v6.2.0.ht
prescriptions_db/cleaned_prescriptions_with_doses_v6.2.0.ht
prescriptions_db/cleaned_prescriptions_splited_to_therapies_v6.2.0.ht
```

The main product of the core preSCRIPT workflow is:

```text
prescriptions_db/cleaned_prescriptions_splited_to_therapies_v6.2.0.ht
```

From this table, users can derive downstream prescription phenotypes such as:

- longest therapy duration,
- median daily dose,
- observed maintenance dose,
- number of therapy episodes,
- dose escalation summaries,
- discontinuation or switching summaries.

Phenotype construction is downstream analysis, not part of the core search-and-annotation workflow.

## Repository structure

```text
prescriptions_processing/
  codes_extractor.py
  drugs_filtering.py
  dose_annotation.py
  data_cleaning.py
  therapies_splitter.py

scripts/
  step_0_prepare_substances_list.ipynb
  step_1_extract_brandnames_and_codes.ipynb
  step_2_fetch_dispensed_ukb_data.ipynb
  step_3_filter_prescriptions.ipynb
  step_4_annotate_dose_and_quantity.ipynb
  step_5_clean_and_impute.ipynb
  step_6_split_into_therapies.ipynb

data/
  input/
    drug_list.csv
    substances.json
    brand_names_refinement_rules.json
  ukb_lkps/
  codes_lkps/
```

## Practical notes

Start with a small drug list when adapting preSCRIPT to a new therapeutic area. Run Steps 0-3 first and inspect the matched prescriptions before investing time in dose annotation and therapy splitting.

If recall is too low, check:

1. `drug_list.csv`,
2. `substances.json`,
3. missing brand names,
4. proper lookup files.

If false positives are too frequent, check:

1. `brand_names_refinement_rules.json`,
2. overly broad alternative names,
3. ambiguous short aliases,
4. brand names shared across products.

The safest way to customize preSCRIPT is to edit the input drug panel and refinement rules, then rerun the workflow from Step 0 or Step 1. The parameters used in `DataCleaning` and `TherapiesSplitter` were selected empirically and can be adjusted for specific research problems. This is especially useful when outlier removal or therapy splitting performs poorly for a new drug panel or dataset.

## Data and licensing notes

### TODO: Licensing related note (Open Gov Lic vs Reserved Rights for dm+d)

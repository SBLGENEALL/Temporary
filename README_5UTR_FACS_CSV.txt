5UTR FACS Dunnett analysis - CSV bundle

FILES
1. 5UTR_FACS_Data_Input.csv
   This is one flat data table. CSV files do not contain multiple sheets.
   Enter values only in Rep1, Rep2, and Rep3.

2. 5UTR_FACS_Dunnett_analysis.R
   Runs one-way ANOVA and two-sided Dunnett comparisons versus Original
   separately for each Day, Selection, and Metric condition.

INPUT RULES
- Keep the header row unchanged.
- Keep the control name exactly as Original.
- Replace TOP1-TOP31 with actual construct names consistently if needed.
- Leave missing measurements blank. Do not use 0 for missing data.
- D1 and D3 are common Pre-selection measurements collected before the split.
- D7 and D11 each contain a w/o Selection arm and a Selection arm.
- Each independently transfected well is one replicate. FACS events are not replicates.
- MFI-like metrics are tested after log2 transformation; plots use the original scale.

RUN DIRECTLY FROM CSV
Keep the CSV and R script in the same server directory, then run:

Rscript 5UTR_FACS_Dunnett_analysis.R --input 5UTR_FACS_Data_Input.csv

CSV-to-TSV conversion is not required. The R script reads CSV directly.

OUTPUT
5UTR_FACS_analysis_results/
  Dunnett_results.csv
  Condition_diagnostics.csv
  Analysis_notes.txt
  figures/

REQUIRED R PACKAGES
dplyr, tidyr, ggplot2, multcomp

The readxl package is not used for CSV input.

If packages are missing, use the company conda mirror:

conda install --solver=libmamba --offline --override-channels \
  -c file:///data/conda_repo/mirror/conda-forge/conda-forge \
  r-dplyr r-tidyr r-ggplot2 r-multcomp

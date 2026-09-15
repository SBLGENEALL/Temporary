5UTR FACS Dunnett analysis - TSV bundle

FILES
1. 5UTR_FACS_Data_Input.tsv
   Enter the three well-level values in Rep1, Rep2, and Rep3.
   The template contains Original and TOP1-TOP31 for these metrics:
   - Live cell (%)
   - GFP positive (%)
   - GFP GeoMean MFI

2. 5UTR_FACS_Dunnett_analysis.R
   Runs one-way ANOVA and two-sided Dunnett comparisons versus Original
   separately for each Day, Selection, and Metric condition.

INPUT RULES
- Keep the control name exactly as Original.
- Replace TOP1-TOP31 with actual construct names consistently if needed.
- Leave missing measurements blank. Do not use 0 to represent missing data.
- D1 is labeled Pre-selection. Later days have Selection and w/o Selection rows.
- MFI-like metrics are tested after log2 transformation, while plots use the original scale.
- Individual FACS events are not replicates. Each independently transfected well is n=1.

RUN
Keep the TSV and R files in the same directory, then run:

Rscript 5UTR_FACS_Dunnett_analysis.R --input 5UTR_FACS_Data_Input.tsv

OUTPUT
5UTR_FACS_analysis_results/
  Dunnett_results.csv
  Condition_diagnostics.csv
  Analysis_notes.txt
  figures/

REQUIRED R PACKAGES
dplyr, tidyr, ggplot2, multcomp

The readxl package is not required when using the supplied TSV input.


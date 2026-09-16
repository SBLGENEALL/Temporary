5UTR FACS Welch-Holm analysis - TSV bundle

FILES
1. 5UTR_FACS_Data_Input.tsv
   Enter the three well-level values in Rep1, Rep2, and Rep3.
   The template contains Original and TOP1-TOP31 for these metrics:
   - Live cell (%)
   - GFP positive (%)
   - GFP GeoMean MFI

2. 5UTR_FACS_Dunnett_analysis.R
   Runs two-sided Welch t-tests versus Original with Holm correction
   separately for each Day, Selection, and Metric condition.
   Dunnett-adjusted results are retained as a secondary comparison.

INPUT RULES
- Keep the control name exactly as Original.
- Replace TOP1-TOP31 with actual construct names consistently if needed.
- Leave missing measurements blank. Do not use 0 to represent missing data.
- D1 and D3 are common Pre-selection measurements collected before the split.
- D7 and D11 each contain a w/o Selection arm and a Selection arm.
- MFI-like metrics are tested after log2 transformation, while plots use the original scale.
- Individual FACS events are not replicates. Each independently transfected well is n=1.

RUN
Keep the TSV and R files in the same directory, then run:

Rscript 5UTR_FACS_Dunnett_analysis.R --input 5UTR_FACS_Data_Input.tsv

OUTPUT
5UTR_FACS_analysis_results/
  Welch_Holm_results.csv
  Dunnett_results.csv (compatibility copy)
  Condition_diagnostics.csv
  Analysis_notes.txt
  figures/

The graph stars and Primary_p_adjusted column use Welch_p_Holm.
Use Welch_Holm_results.csv as the primary result table.

REQUIRED R PACKAGES
dplyr, tidyr, ggplot2, multcomp

The readxl package is not required when using the supplied TSV input.

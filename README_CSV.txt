CHO 5'UTR Transposase FACS CSV bundle

Files
1. Plate_Map.csv
   - Enter Construct, Replicate, DNA_Prep, DNA_Batch once for wells A1-H12.
   - Plot order follows the first appearance of each Construct in Plate_Map; TOP numbers are not sorted.

2. FACS_Data.csv
   - Paste FlowJo values into the matching Day and Well rows.
   - Condition values: Pre_split, GS, No_selection
   - Enter 50% as 50, not 0.50.

3. plot_facs_longitudinal.R
   - The script accepts CSV or TSV pairs and prefers CSV when both exist.
   - Run: Rscript plot_facs_longitudinal.R
   - Results are written to FACS_plot_results_v3_prism/.
   - Main plots:
       01_MFI_barplots_by_day/ (one mean +/- SD bar plot per day and condition)
       02_GFP_positive_barplots_by_day/ (one mean +/- SD bar plot per day and condition)
       03_MFI_lineplot_top_candidates.png
       04_GFP_positive_lineplot_top_candidates.png
       05_MFI_lineplot_all_candidates_small_multiples.png

Recommended workflow
- Transfer this ZIP to the internal Linux server without opening its contents in Excel first.
- Extract on Linux:
    unzip CHO_5UTR_FACS_CSV_Bundle.zip
    cd CHO_5UTR_FACS_CSV_Bundle
- Edit the CSV files using LibreOffice Calc and save as UTF-8, comma-delimited CSV.
- If editing as plain text, keep commas between columns; do not paste tab-delimited rows directly.
- Keep the exact filenames shown above.
- Run:
    Rscript plot_facs_longitudinal.R

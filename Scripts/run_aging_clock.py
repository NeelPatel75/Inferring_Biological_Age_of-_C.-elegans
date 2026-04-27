# All necesary modules 
import os
import pandas as pd
import glob
import sys
import numpy as np
from scipy.stats import ttest_ind

# This py file will run the entire pipeline using the kallisto outputs (Phase 2-4)
# The results can be found in the Results folder under the names: 
# GSE65765_CPM.csv for phase 2 
# Final_Biological_Age_Results.csv for phase 3 
# Final_Statistical_Analysis.txt for phase 4 
# Make sure the results from run_kallisto.py (Phase 1) are located in the Results folder
# under the name kallisto_outputs before running this file.
# If there is an error while running this pipline please refer to the print statments in
# the terminal to see what the problem is. 
#______________________________________________________________________

#-- 1. Set up File Paths --
# The file path set up is generic therfore the user can simply run this py file for the pipline to run 
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
# PROJECT_ROOT is the parent folder (C.-elegans/)
PROJECT_ROOT = os.path.dirname(SCRIPT_DIR)
# Path to all of the data 
DATA_DIR = os.path.join(PROJECT_ROOT, "Data")
# Path so that all of the results will end up in the correct folder 
RESULTS_DIR = os.path.join(PROJECT_ROOT, "Results")
# Path to the kallisto results 
KALLISTO_DIR = os.path.join(RESULTS_DIR, "kallisto_outputs")
# Ensure all folders exits 
os.makedirs(RESULTS_DIR, exist_ok=True)
# Define file paths
# File path for results of phase 2 
CPM_FILE = os.path.join(RESULTS_DIR, "GSE65765_CPM.csv")
# File path for  results of phase 3 
OUTPUT_FILE = os.path.join(RESULTS_DIR, "Final_Biological_Age_Results.csv")
# File path for results of phase 4 
STATS_REPORT_TXT = os.path.join(RESULTS_DIR, "Final_Statistical_Analysis.txt")
# File path to the predictor genes used in phase 2/3 
PREDICTOR_FILE = os.path.join(DATA_DIR, "Predictor_Genes.csv")

# Add SCRIPT_DIR to path for module imports
if SCRIPT_DIR not in sys.path:
    sys.path.append(SCRIPT_DIR)
# Call in the scripts from the BiT Age Tool 
try:
    from Biological_Age_Predictor import predict 
    from biological_age_correction import calculate_Bio_Age_correction
except ModuleNotFoundError as e:
    sys.exit(f"ERROR: {e}. Ensure helper scripts are in {SCRIPT_DIR}")

# --- 2. validate predictor files  ---
# Check if the Predictor_Genes.csv file actually exists at the defined path
# This prevents the script from crashing with a confusing pandas error later
if not os.path.exists(PREDICTOR_FILE):
    # Tells the user the file is missing 
    sys.exit(f"ERROR: Predictor file not found at {PREDICTOR_FILE}")
# Load the predictor file into a DataFrame
# index_col=0 tells pandas to use the first column (the Gene IDs) as the row labels
predictor_df = pd.read_csv(PREDICTOR_FILE, index_col=0)
# Extract the gene IDs 
required_genes = predictor_df.index

# --- 3. Aggregation Phase (Phase 2 code) ---
# Tells the user that phase 2 is starting 
print(f"\nStep 1: Aggregating samples from {KALLISTO_DIR}...")
# Search pattern to find all files ending in '_final.csv' within the kallisto_outputs 
search_path = os.path.join(KALLISTO_DIR, "**", "*_final.csv")
file_list = glob.glob(search_path, recursive=True)
# If no files are found, exit the script to avoid errors in the merge step
# Tell the user the files cannot be found 
if not file_list:
    sys.exit(f"ERROR: No '_final.csv' files found in {KALLISTO_DIR}")
# This list will hold the processed data for each individual sample
all_samples = []
# Loop through every identified file to process and clean the transcript data
for file in file_list:
    # Extract the SRR ID from the filename
    sample_name = os.path.basename(file).replace("_final.csv", "")
    # Load the specific sample file 
    df_sample = pd.read_csv(file)
    # Skip the file if it doesn't contain the necessary Gene ID column
    if 'WBGene_ID' not in df_sample.columns: continue
    # Collapse transcript-level data into gene-level data
    # 1. Remove rows where the Gene ID is missing (NaN)
    # 2. Group by the Gene ID and sum the TPM values of all related transcripts
    gene_level_df = df_sample.dropna(subset=['WBGene_ID']).groupby('WBGene_ID')['tpm'].sum().reset_index()
    # Label the columns: 'Gene_ID' and the specific SRR sample name
    gene_level_df.columns = ['Gene_ID', sample_name]
    # Set Gene_ID as the index to allow for easy alignment/merging with other samples
    all_samples.append(gene_level_df.set_index('Gene_ID'))
# Merge all individual sample DataFrames into one large 'Master' table
# Aligning on 'axis=1' ensures samples are columns and genes are rows
# .reindex ensures we only keep the genes required by the predictor
# .fillna(0) ensures that if a gene was missing in a sample, it is treated as 0 expression
master_df = pd.concat(all_samples, axis=1).reindex(required_genes).fillna(0)
# Convert TPM (Transcripts Per Million) to CPM (Counts Per Million)
# Divide each value by the sample's total expression and multiply by 1,000,000
master_df = master_df.divide(master_df.sum(axis=0), axis=1) * 1e6
master_df.to_csv(CPM_FILE, sep='\t')

# --- 4. Run the BiT Age Tool (Phase 3) ---
# Call the predict function from Biological_Age_Predictor.py
# This uses the machine learning model to estimate age based on the CPM matrix
print("\nStep 2 & 3: Running Prediction and Correction...")
raw_results = predict(CPM_FILE, PREDICTOR_FILE)
# Standardize the output column name to 'Bio_Age' for compatibility with the correction script
raw_results = raw_results.rename(columns={'Predicted_Biological_Age': 'Bio_Age'})
# Apply the statistical survival correction to the raw predicted ages
# This adjusts for biases in the biological clock models
final_df = calculate_Bio_Age_correction(raw_results)
# Save the final corrected biological age results to a CSV file
final_df.to_csv(OUTPUT_FILE)
# --- 5. Test for statistical significance using T-test (Phase 4) ---
print("\nStep 4: Running Detailed Statistical Analysis...")

# Moves SRR IDs from index to a column named Sample to avoid KeyError
if 'Sample' not in final_df.columns:
    final_df = final_df.reset_index()
    final_df = final_df.rename(columns={final_df.columns[0]: 'Sample'})

# Define the gorupings for the SRR files 
wt_embryo = ["SRR28479540", "SRR28479539", "SRR28479538", "SRR28479525", "SRR28479526", "SRR28479527"]
wt_l1 = ["SRR28479534", "SRR28479533", "SRR28479532", "SRR28479531"]
xol1_embryo = ["SRR28479537", "SRR28479536", "SRR28479535"]
xol1_l1 = ["SRR28479530", "SRR28479529", "SRR28479528"]
# Create aggregated lists for all of the necessary comparisons
all_embryo = wt_embryo + xol1_embryo
all_l1 = wt_l1 + xol1_l1
wt_all = wt_embryo + wt_l1
xol1_all = xol1_embryo + xol1_l1

# Finds and extracts the biological age values for a specific list of SRR IDs.
def get_group_ages(df, srr_list, sample_col, age_col):
    # Uses .startswith() to ensure matches even if Sample IDs have suffixes
    mask = df[sample_col].apply(lambda x: any(str(x).startswith(srr) for srr in srr_list))
    return df[mask][age_col].dropna()
# Performs Welch's T-Test
def run_t_test(group1, group2, label1, label2, file=None):
    # A T-test requires at least two samples in each group to calculate variance
    # If not satisfied tells user 
    if len(group1) < 2 or len(group2) < 2:
        msg = f"\n[SKIP] {label1} vs {label2}: Insufficient samples (n1={len(group1)}, n2={len(group2)})"
        print(msg)
        if file: file.write(msg + "\n")
        return
    # Calculate t-statistic and p-value
    t_stat, p_val = ttest_ind(group1, group2, equal_var=False)
    # Format the results into a readable list of strings
    res = [
        f"\nCOMPARISON: {label1} vs {label2}",
        f"  Sample sizes: n({label1})={len(group1)}, n({label2})={len(group2)}",
        f"  Mean Ages:    {label1}={group1.mean():.4f}, {label2}={group2.mean():.4f}",
        f"  t={t_stat:.4f}, p={p_val:.6f}",
        f"  SIGNIFICANT:  {'YES' if p_val < 0.05 else 'NO'}"
    ]
    # If the result is significant (p < 0.05), provide a biological interpretation
    if p_val < 0.05:
        direction = "OLDER" if group1.mean() > group2.mean() else "YOUNGER"
        res.append(f"  INTERPRETATION: {label1} is biologically {direction} than {label2}")
    # Output the results to both the terminal and the report file
    for line in res:
        print(line)
        if file: file.write(line + "\n")

# EXECUTION OF T-TESTS
sample_col, age_col = "Sample", "Corrected_Biological_Age"

# Separate out the specific age data for each subset defined above
ages_wt_embryo = get_group_ages(final_df, wt_embryo, sample_col, age_col)
ages_wt_l1 = get_group_ages(final_df, wt_l1, sample_col, age_col)
ages_xol1_embryo = get_group_ages(final_df, xol1_embryo, sample_col, age_col)
ages_xol1_l1 = get_group_ages(final_df, xol1_l1, sample_col, age_col)

# Separate out the combined age data for the different comparisons 
ages_all_embryo = get_group_ages(final_df, all_embryo, sample_col, age_col)
ages_all_l1 = get_group_ages(final_df, all_l1, sample_col, age_col)
ages_wt_all = get_group_ages(final_df, wt_all, sample_col, age_col)
ages_xol1_all = get_group_ages(final_df, xol1_all, sample_col, age_col)

with open(STATS_REPORT_TXT, "w") as f:
    f.write("C. ELEGANS BIOLOGICAL AGE STATISTICAL REPORT\n" + "="*60 + "\n")

    # Control comparison 
    # Checks if the clock correctly sees L1 as older than Embryos
    # If this is false it suggests that there is an error in the pipeline
    f.write("\n--- SECTION 1: CLOCK VALIDATION (WT ONLY) ---\n")
    run_t_test(ages_wt_l1, ages_wt_embryo, "WT_L1", "WT_Embryo", file=f)

    # General State Analysis 
    # Compares the average age of all larvae vs all embryos.
    f.write("\n--- SECTION 2: GENERAL DEVELOPMENTAL STAGES (ALL SAMPLES) ---\n")
    run_t_test(ages_all_l1, ages_all_embryo, "All_L1_Larvae", "All_Embryos", file=f)

    # Compare xol-1 mutants with WT within the same stage 
    f.write("\n--- SECTION 3: GENOTYPE IMPACT WITHIN STAGES ---\n")
    # Early Embryo
    run_t_test(ages_xol1_embryo, ages_wt_embryo, "xol1_Embryo", "WT_Embryo", file=f)
    # L1
    run_t_test(ages_xol1_l1, ages_wt_l1, "xol1_L1", "WT_L1", file=f)
    # General test of xol-1 mutants vs WT 
    f.write("\n--- SECTION 4: OVERALL IMPACT (All Stages Combined) ---\n")
    run_t_test(ages_xol1_all, ages_wt_all, "xol1_All", "WT_All", file=f)
# Let the user know the pipling worked correctly 
print(f"\n--- Workflow Complete. Report found at: {STATS_REPORT_TXT} ---")
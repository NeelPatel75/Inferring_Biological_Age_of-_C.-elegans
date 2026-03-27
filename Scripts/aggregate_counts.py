import os
import pandas as pd

# Path to the pre-built Kallisto index file
# This index is created from the C. elegans reference transcriptome
BASE_DIR = "/home/npatel79/C_Elegans_COMP"
INPUT_DIR = os.path.join(BASE_DIR, "Results/Kallisto_Outputs")
# Even though it's called .csv, we save it with tabs to match your tool's default
OUTPUT_FILE = os.path.join(BASE_DIR, "GSE65765_CPM.csv")

# Generate the list of 16 samples
samples = [f"SRR284795{i}" for i in range(25, 41)]
all_data = []

print("Aggregating Kallisto outputs for BiT Tool compatibility...")

for srr in samples:
    file_path = os.path.join(INPUT_DIR, srr, "abundance.tsv")
    
    if os.path.exists(file_path):
        # Read the raw counts
        df = pd.read_csv(file_path, sep='\t', usecols=['target_id', 'est_counts'])
        
        # --- CRITICAL STEP: CLEAN IDS ---
        # Remove transcript version numbers (e.g., WBGene00000001.1 -> WBGene00000001)
        # Your tool needs the 'WBG' prefix to find the genes in the regex filter.
        df['target_id'] = df['target_id'].str.split('.').str[0]
        
        # Rename column to the Sample ID
        df = df.rename(columns={'est_counts': srr})
        df.set_index('target_id', inplace=True)
        
        all_data.append(df)
        print(f"Processed: {srr}")
    else:
        print(f"Warning: Missing file for {srr}")

# 2. Combine all 16 samples
master_df = pd.concat(all_data, axis=1)

# 3. CPM Normalization
# (Counts / Total Reads) * 1,000,000
for srr in samples:
    total_reads = master_df[srr].sum()
    if total_reads > 0:
        master_df[srr] = (master_df[srr] / total_reads) * 1e6

# 4. Final Formatting for the 'predict' function
# Name the index column exactly 'Gene_ID'
master_df.index.name = 'Gene_ID'

# 5. Save as TAB-SEPARATED (to match the 'sep' in your tool code)
master_df.to_csv(OUTPUT_FILE, sep='\t')

print(f"\nPhase 2 Complete! File saved: {OUTPUT_FILE}")
print(f"Shape: {master_df.shape} (Genes x Samples)")

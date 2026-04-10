import os
import pandas as pd
import glob

# --- Configuration ---
BASE_DIR = "/home/npatel79/Inferring_Biological_Age_of-_C.-elegans"
INPUT_DIR = os.path.join(BASE_DIR, "kallisto_outputs")
PREDICTOR_FILE = os.path.join(BASE_DIR, "Predictor_Genes.csv") # Path to your coefficients file
FINAL_OUTPUT_PATH = os.path.join(BASE_DIR, "GSE65765_CPM.csv")

# 1. Load the required genes from the predictor file to prevent KeyErrors
print("Loading required predictor genes...")
predictor_df = pd.read_csv(PREDICTOR_FILE, index_col=0)
required_genes = predictor_df.index

# 2. Find all sample-specific final.csv files created in Phase 1
search_path = os.path.join(INPUT_DIR, "**", "*_final.csv")
file_list = glob.glob(search_path, recursive=True)

if not file_list:
    print("No '_final.csv' files found. Please ensure Phase 1 finished successfully.")
else:
    all_samples = []

    for file in file_list:
        sample_name = os.path.basename(file).replace("_final.csv", "")
        print(f"Processing: {sample_name}")
        
        df = pd.read_csv(file)
        
        # Filter and aggregate transcript TPMs to Gene level
        gene_level_df = df.dropna(subset=['WBGene_ID']).groupby('WBGene_ID')['tpm'].sum().reset_index()
        
        gene_level_df.columns = ['Gene_ID', sample_name]
        gene_level_df = gene_level_df.set_index('Gene_ID')
        
        all_samples.append(gene_level_df)

    # 3. Merge all samples
    print("Merging samples and aligning with predictor genes...")
    master_df = pd.concat(all_samples, axis=1)

    # 4. FIX: Reindex ensures all 576 required genes exist in the table.
    # If a gene was missing from your data, it is added here with a value of 0.
    master_df = master_df.reindex(required_genes).fillna(0)
    
    # 5. Save as Tab-Separated Values
    master_df.to_csv(FINAL_OUTPUT_PATH, sep='\t')
    print(f"Done! Master file saved to: {FINAL_OUTPUT_PATH}")
    print(f"The file now contains {len(master_df)} genes, matching the predictor's requirements.")
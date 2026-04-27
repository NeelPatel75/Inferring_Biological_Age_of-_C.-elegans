# All necesary modules 
import os
import subprocess
import pandas as pd
import glob
import sys

# This script will create an abundance.tsv for each SRR which will contain 
# all of the abundances for each gene. 
# This Script will also create a SRRXXXX_final.csv for each SRR which will match each 
# target_id from the abundance.tsv to its WormBase ID for phase 2-4. 
# Note: in order to use this script you will first need to download the necessary indexes 
# from https://parasite.wormbase.org/Caenorhabditis_elegans_prjna13758/Info/Index/ 
# in order to run this script. 
# This git hub contains a transcript_map.csv that analyses the files found on that 
# website so that the user does not have to download them. 
# Instructions on how to run the script so can be found in the README of the Git Hub. 
#____________________________________________________________________________

# --- 1. GENERIC DIRECTORY SETUP ---
# --- 1. GENERIC DIRECTORY SETUP ---
# Get the absolute path of the directory where THIS script is located (e.g., /.../Scripts/)
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))

# Go up one level to the project root (e.g., /.../Inferring_Biological_Age_of-_C.-elegans/)
BASE_DIR = os.path.dirname(SCRIPT_DIR)

# Point this to where your FASTQs and FASTA are currently located
DATA_DIR = os.path.join(BASE_DIR, "Data", "sample_test_data")
OUTPUT_DIR = os.path.join(BASE_DIR, "Results", "kallisto_outputs")
os.makedirs(OUTPUT_DIR, exist_ok=True)

# --- 2. UPDATED FILE PATHS ---
# Using the exact name of the healthy file we just verified
TRANSCRIPTOME_FA = os.path.join(DATA_DIR, "caenorhabditis_elegans.PRJNA13758.WBPS19.mRNA_transcripts.fa")
INDEX_FILE = os.path.join(BASE_DIR, "celegans_index.idx")
MAPPING_FILE = os.path.join(BASE_DIR, "Data", "transcript_map.csv")

# --- 3. INDEX BUILDING ---
if not os.path.exists(INDEX_FILE):
    if os.path.exists(TRANSCRIPTOME_FA):
        print(f"Building index...")
        subprocess.run(["kallisto", "index", "-i", INDEX_FILE, TRANSCRIPTOME_FA], check=True)
    else:
        sys.exit(f"ERROR: Index and FASTA missing.")

# --- 4. KALLISTO QUANTIFICATION ---
# Using DATA_DIR because that is where we defined the path earlier
fastq_files = sorted(glob.glob(os.path.join(DATA_DIR, "*.fastq")))
samples = {}
for fq in fastq_files:
    name = os.path.basename(fq)
    if "_1.fastq" in name:
        base = name.replace("_1.fastq", "")
        # Again, using DATA_DIR here
        mate = os.path.join(DATA_DIR, base + "_2.fastq")
        if os.path.exists(mate):
            samples[base] = (fq, mate)

for sample, (fq1, fq2) in samples.items():
    out_dir = os.path.join(OUTPUT_DIR, sample)
    if not os.path.exists(out_dir):
        print(f"Processing {sample}...")
        os.makedirs(out_dir, exist_ok=True)
        subprocess.run(["kallisto", "quant", "-i", INDEX_FILE, "-o", out_dir, "-b", "100", fq1, fq2], check=True)

# --- 5. MAPPING USING THE LIGHTWEIGHT CSV ---
# Using MAPPING_FILE as defined in Section 2
if not os.path.exists(MAPPING_FILE):
    sys.exit(f"ERROR: {MAPPING_FILE} not found. Ensure it is in Data/transcript_map.csv")

print(f"Loading transcript-to-gene map from {os.path.basename(MAPPING_FILE)}...")
map_df = pd.read_csv(MAPPING_FILE)

for sample in samples.keys():
    abundance_path = os.path.join(OUTPUT_DIR, sample, "abundance.tsv")
    if os.path.exists(abundance_path):
        df = pd.read_csv(abundance_path, sep='\t')
        
        # Clean target IDs to match the Map (strips "transcript:" prefix if present)
        df['target_id'] = df['target_id'].astype(str).str.split(":").str[-1]
        
        # Merge with our pre-created map
        final_df = df.merge(map_df, on="target_id", how="left")
        
        output_file = os.path.join(OUTPUT_DIR, sample, f"{sample}_final.csv")
        final_df[['target_id', 'WBGene_ID', 'tpm']].to_csv(output_file, index=False)
        print(f"Sample {sample}: Final CSV created.")

print("\nPhase 1 execution finished.")
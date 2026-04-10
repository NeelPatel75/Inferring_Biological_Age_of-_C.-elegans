import os
import subprocess
import pandas as pd
import glob

# This script will create an abundance.tsv for each SRR which will contain 
# all of the abundances for each gene. 
# This Script will also create a SRRXXXX_final.csv for each SRR which will match each 
# target_id from the abundance.tsv to its WormBase ID for phase 2-4. 
# Note: in order to use this script you will first need to download teh necessary indexes 
# from https://parasite.wormbase.org/Caenorhabditis_elegans_prjna13758/Info/Index/ 
# in order to run this script. 
# Instructions on how to do so can be found in the README of the Git Hub. 
#____________________________________________________________________________

# Start of Script 
# Define the root directory for the C. elegans aging project
BASE_DIR = "/home/npatel79/Inferring_Biological_Age_of-_C.-elegans"
DATA_DIR = os.path.join(BASE_DIR, "sample_test_data")
OUTPUT_DIR = os.path.join(BASE_DIR, "kallisto_outputs")
# Kallisto uses a binary 'index' of the transcriptome to map reads
INDEX_FILE = os.path.join(BASE_DIR, "celegans_index.idx")
# Use glob to find reference files
fasta_search = glob.glob(os.path.join(DATA_DIR, "*mRNA_transcripts.fa"))
gff_search = glob.glob(os.path.join(DATA_DIR, "*annotations.gff3"))
# Select the first match for the FASTA (Transcriptome) and GFF (Annotations)
TRANSCRIPTOME_FA = fasta_search[0] if fasta_search else None
GFF_FILE = gff_search[0] if gff_search else None
# Make sure the output file exists 
os.makedirs(OUTPUT_DIR, exist_ok=True)

# Turns the FASTA file into a searchable index.
if TRANSCRIPTOME_FA and not os.path.exists(INDEX_FILE):
    print(f"Building index from {os.path.basename(TRANSCRIPTOME_FA)}...")
    # kallisto index creates the .idx file required for the quant step
    subprocess.run(["kallisto", "index", "-i", INDEX_FILE, TRANSCRIPTOME_FA])

# Looks for sample through the suffix 
fastq_files = sorted(glob.glob(os.path.join(DATA_DIR, "*.fastq")))
samples = {}
for fq in fastq_files:
    name = os.path.basename(fq)
    if "_1.fastq" in name:
        base = name.replace("_1.fastq", "")
        mate = os.path.join(DATA_DIR, base + "_2.fastq")
        # Store as a tuple (Forward, Reverse) if the second file exists
        if os.path.exists(mate):
            samples[base] = (fq, mate)

# Maps the raw reads to the transcripts and calculates abundance (TPM/Counts)
for sample, (fq1, fq2) in samples.items():
    out_dir = os.path.join(OUTPUT_DIR, sample)
    # Performs 100 bootstrap samples for better statistical reliability
    if not os.path.exists(out_dir):
        print(f"Processing {sample}...")
        subprocess.run(["kallisto", "quant", "-i", INDEX_FILE, "-o", out_dir, "-b", "100", fq1, fq2])
# Uses the GFF3 file to build a map with the Kallisto outputs Transcript IDs (isoforms) 
# As a result it gives the WormBase ID in the for the out puts as needed by the Bit Age tool. 
mapping = []
if GFF_FILE:
    print(f"--- Processing GFF3: {os.path.basename(GFF_FILE)} ---")
    with open(GFF_FILE) as f:
        for line in f:
            # Filter for mRNA/transcript rows which contain both ID and Parent Gene
            if "\tmRNA\t" in line or "\ttranscript\t" in line:
                attrs = line.split("\t")[8]
                if "ID=" in attrs and "Parent=" in attrs:
                    # Split attributes by semicolon
                    parts = {item.split("=")[0]: item.split("=")[1] for item in attrs.split(";") if "=" in item}
                    
                    # Extract and clean IDs for proper formatting 
                    t_id = parts.get("ID", "").split(":")[-1]
                    g_id = parts.get("Parent", "").split(":")[-1]
                    
                    if t_id and g_id:
                        mapping.append({"target_id": t_id, "WBGene_ID": g_id})
# Convert the list to a DataFrame for easy merging
map_df = pd.DataFrame(mapping).drop_duplicates()
print(f"--- Successfully created a map for {len(map_df)} genes ---")

# Combines the abundance data with the WBGene IDs for a final clean output
for sample in samples.keys():
    abundance_path = os.path.join(OUTPUT_DIR, sample, "abundance.tsv")
    if os.path.exists(abundance_path):
        # Load the raw Kallisto output
        df = pd.read_csv(abundance_path, sep='\t')
        
        # Clean transcript IDs in the TSV to match the cleaned IDs from the GFF
        df['target_id'] = df['target_id'].astype(str).str.split(":").str[-1]
        
        # Perform a left merge to attach WBGene_ID to each transcript row
        final_df = df.merge(map_df, on="target_id", how="left")
        
        # Validation print to see if mapping was successful for this sample
        num_matched = final_df['WBGene_ID'].notna().sum()
        print(f"Sample {sample}: Matched {num_matched} / {len(df)} transcripts.")
        
        # Export 'target_id', 'WBGene_ID', and 'tpm' (required for Phase 2 aggregation)
        output_file = os.path.join(OUTPUT_DIR, sample, f"{sample}_final.csv")
        final_df[['target_id', 'WBGene_ID', 'tpm']].to_csv(output_file, index=False)

print("\nPhase 1 execution finished.")
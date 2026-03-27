import os

# Path to the pre-built Kallisto index file
# This index is created from the C. elegans reference transcriptome
INDEX = "Data/elegans_index.idx"
# Directory containing all downloaded RNA-seq FASTQ files
# Each sample should have paired-end reads: *_1.fastq and *_2.fastq
DATA_DIR = "Transcriptome_data"
# Output directory where Kallisto results will be stored
OUT_DIR = "Results/Kallisto_Outputs"

# Create the output directory if it does not already exist
os.makedirs(OUT_DIR, exist_ok=True)
# Loop through the SRR IDs from SRR28479525 to SRR28479540
for i in range(25, 41):
     # Construct the full SRR accession number
    srr = f"SRR284795{i}"
    print(f"\n>>> Mapping {srr} <<<")
    
    # Input files are inside Transcriptome_data
    # Paired-end read files for each sample
    r1 = f"{DATA_DIR}/{srr}_1.fastq"
    r2 = f"{DATA_DIR}/{srr}_2.fastq"
     # Create a separate folder for each sample's results
    sample_out = f"{OUT_DIR}/{srr}"
    os.makedirs(sample_out, exist_ok=True)
    
    # Run the kallisto quant function to quantify the data 
    # Kallisto performs pseudo-alignment and outputs:
    # abundance.tsv (gene expression values) (I will rename soon to be more specific)
    # run_info.json (run metadata)
    os.system(f"kallisto quant -i {INDEX} -o {sample_out} {r1} {r2}")

print("\n--- Done! Check Results/Kallisto_Outputs for your results ---")
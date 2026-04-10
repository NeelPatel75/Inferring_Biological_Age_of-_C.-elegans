import os
# 1. Setup Paths
INDEX = "/home/npatel79/Inferring_Biological_Age_of-_C.-elegans/Data/elegans_index.idx"
DATA_DIR = "/home/npatel79/Inferring_Biological_Age_of-_C.-elegans/sample_test_data"
OUT_DIR = "/home/npatel79/Inferring_Biological_Age_of-_C.-elegans/test_results/Kallisto_samples"

os.makedirs(OUT_DIR, exist_ok=True)

# 2. UPDATE THIS RANGE: Only use 25 and 26 for the test
# range(25, 27) will run for SRR28479525 and SRR28479526
for i in range(25, 27): 
    srr = f"SRR284795{i}"
    print(f"\n>>> Mapping {srr} <<<")
    
    r1 = f"{DATA_DIR}/{srr}_1.fastq"
    r2 = f"{DATA_DIR}/{srr}_2.fastq"
    
    # 3. SAFETY CHECK: Only run if the files actually exist
    if os.path.exists(r1) and os.path.exists(r2):
        sample_out = f"{OUT_DIR}/{srr}"
        os.makedirs(sample_out, exist_ok=True)
        
        # Execute Kallisto
        os.system(f"kallisto quant -i {INDEX} -o {sample_out} {r1} {r2}")
    else:
        print(f"Skipping {srr}: Fastq files not found in {DATA_DIR}")

print(f"\n--- Done! Check {OUT_DIR} for your results ---")
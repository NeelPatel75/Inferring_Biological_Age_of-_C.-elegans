# Inferring Biological Age of *C. elegans* via Transcriptomic Aging Clocks

## Overview

This project implements a bioinformatics pipeline to estimate the **biological age** of *Caenorhabditis elegans* using RNA-seq data.

The workflow integrates:

- **Kallisto** for transcript quantification
- A **biological aging clock model** for age prediction
- **Welch’s t-tests** for statistical analysis

This pipeline allows users to study how genetic mutations, such as *xol-1*, impact developmental timing and biological aging.

---

## Project Structure

| File | Role |
|------|------|
| `run_kallisto.py` | Runs Phase 1 and creates the Kallisto output folder |
| `run_aging_clock.py` | Runs Phases 2–4 using the Kallisto outputs |
| `Biological_Age_Predictor.py` | Aging clock model (called automatically by run_aging_clock.py) |
| `biological_age_correction.py` | Correction step (called automatically by run_aging_clock.py) |
| `mapping_file.py` | Used to map WBgne_ID with target_ID (already ran to set up for phase 1, user does not need to run this) |



## Pipeline Overview

### Phase 1: Read Quantification with Kallisto

- Builds the transcriptome index
- Processes paired-end FASTQ files
- Generates abundance values (TPM)
- Maps transcripts to WormBase Gene IDs

### Phase 2: Data Aggregation and Normalization

- Converts transcript-level data to gene-level data
- Aggregates TPM values across samples
- Converts TPM values into CPM values

### Phase 3: Biological Age Prediction

- Uses the predictor gene set
- Applies the biological aging clock model
- Outputs predicted and corrected biological age values

### Phase 4: Statistical Analysis

- Performs **Welch’s t-tests**
- Compares:
  - WT Embryo vs WT L1
  - xol-1 Embryo vs WT Embryo
  - xol-1 L1 vs WT L1
  - All xol-1 vs All WT
- Outputs p-values, group means, and interpretation

---

## Dependencies

### Software/Modules used 

- Python 3.8+
- Kallisto
- Pandas
- NumPy
- Spicy Stats
- glob
- sys
- subprocess
- os 

### How to Download Python Packages 

```bash
pip install pandas numpy scipy statsmodels pingouin
```

---

## Installation

### 1. Clone the Repository

```bash
git clone https://github.com/your-username/Inferring_Biological_Age_of_C_elegans.git
cd Inferring_Biological_Age_of_C_elegans
```

### 2. Install Kallisto

```bash
conda install -c bioconda kallisto
kallisto version 

```
---

## How to Run the Pipeline

### Step 1: Run Kallisto

```bash
python3 Scripts/run_kallisto.py
```

### Step 2: Run Aging Clock Pipeline

```bash
python3 Scripts/run_pipeline.py
```

---

## Output Files

- **kallisto_outputs**.csv → Abundances in TPM + WBGene map 
- **GSE65765_CPM.csv** → gene expression matrix in CPM 
- **Final_Biological_Age_Results.csv** → predicted biological age
- **Final_Statistical_Analysis.txt** → statistical results

##### Note: Look at results folder for all of these files 

---

## Sample Test Data

Reduced-size RNA-seq data from NCBI SRA for fast testing.
Note: Due to this dataset being reduced, do not take any of the results of the sample test data as an accurate representation of biological age/ analysis. This sample data is simply to be used to run the pipeline and analyze how it works. 

---

## Notes

- `run_kallisto.py` works with new datasets and is phase 1 
- `run_again_clock.py` runs the full pipeline phase 2-4 
- Files are automatically detected
- Ensure correct file paths

---

## Troubleshooting

- FASTQ not detected → check folder + naming
- Missing predictor → `Data/Predictor_Genes.csv`
- Kallisto error → run `kallisto version`
- Look at print statements in the terminal to check where the error occurred.  

---

## Future Improvements

- Add visualizations
- Support larger datasets
- Biological data improvement: look at the Adult C. Elegans transcriptome

---

## Citations

- https://pachterlab.github.io/kallisto/manual
- https://docs.scipy.org/doc/scipy/reference/stats.html
- https://docs.python.org/3/library/glob.html
- https://www.alliancegenome.org/members/wormbase
- https://www.dataquest.io/blog/python-subprocess/
- https://www.geeksforgeeks.org/python/python-programming-language-tutorial/

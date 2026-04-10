import pandas as pd
# Import the functions from your existing files
# Ensure these filenames match your .py files exactly (minus the .py)
from Biological_Age_Predictor import predict 
from biological_age_correction import calculate_Bio_Age_correction

# --- SETTINGS ---
cpm_file = '/home/npatel79/Inferring_Biological_Age_of-_C.-elegans/GSE65765_CPM.csv'
genes_file = '/home/npatel79/Inferring_Biological_Age_of-_C.-elegans/Predictor_Genes.csv'
output_file = 'Final_Biological_Age_Results.csv'

print("--- Starting Biological Age Prediction Workflow ---")

# 1. Run the prediction (Part 1)
# This uses your first script to get the initial age estimate
raw_results = predict(cpm_file, genes_file)

# 2. Rename the column
# We must change 'Predicted_Biological_Age' to 'Bio_Age' 
# so the second script can find the data it needs.
raw_results = raw_results.rename(columns={'Predicted_Biological_Age': 'Bio_Age'})

# 3. Run the correction (Part 2)
# This uses your second script to apply the statistical survival adjustment
final_df = calculate_Bio_Age_correction(raw_results)

# 4. Save the results
final_df.to_csv(output_file)

print(f"Success! The final corrected ages are saved in: {output_file}")
print("Preview of results:")
print(final_df.head())
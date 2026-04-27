import pandas as pd

# Path to your giant GFF3
GFF_FILE = "/home/npatel79/Inferring_Biological_Age_of-_C.-elegans/Data/sample_test_data/caenorhabditis_elegans.PRJNA13758.WBPS19.annotations.gff3"

mapping = []
print("Scanning GFF3 to create a lightweight map...")
with open(GFF_FILE) as f:
    for line in f:
        if "\tmRNA\t" in line or "\ttranscript\t" in line:
            attrs = line.split("\t")[8]
            if "ID=" in attrs and "Parent=" in attrs:
                parts = {item.split("=")[0]: item.split("=")[1] for item in attrs.split(";") if "=" in item}
                t_id = parts.get("ID", "").split(":")[-1]
                g_id = parts.get("Parent", "").split(":")[-1]
                if t_id and g_id:
                    mapping.append({"target_id": t_id, "WBGene_ID": g_id})

map_df = pd.DataFrame(mapping).drop_duplicates()
# Save this small file to your Data folder
map_df.to_csv("Data/transcript_map.csv", index=False)
print("Done! You can now upload 'Data/transcript_map.csv' to GitHub.")
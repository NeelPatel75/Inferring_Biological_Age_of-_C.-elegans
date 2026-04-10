import pandas as pd

# 1. Load your current file
df = pd.read_csv('GSE65765_CPM.csv', sep='\t', index_col=0)

# 2. Check if the index contains WBGene. If NOT, flip it.
if not str(df.index[0]).startswith('WBG'):
    print("Detected transposed data. Flipping now...")
    df = df.T 

# 3. Ensure the index name is 'Gene_ID'
df.index.name = 'Gene_ID'

# 4. Clean the index (Remove .1 version numbers and spaces)
df.index = df.index.str.split('.').str[0].str.strip()

# 5. Save it back
df.to_csv('GSE65765_CPM.csv', sep='\t')

print("Fixed! Your Gene IDs are now in the correct column.")
print("New Index Check:", df.index[:3].tolist())
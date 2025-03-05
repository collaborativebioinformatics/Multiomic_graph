import pandas as pd
import re

### USER INPUTS ###
intersect_file = "test/intersected_results.bed"
pqtl_combined_file = "test/pqtl_combined.tsv"
gene_list_file = "test/test_gene_list.csv"
output_file = "test/final_filtered_pqtl.tsv"

print("🔹 Loading datasets...")

# Load `intersected_results.bed` **with strict TAB separation**
intersected_df = pd.read_csv(intersect_file, sep="\t", header=None, dtype=str)

# Debugging: Print number of columns detected
print(f"🔹 Detected {intersected_df.shape[1]} columns in `intersected_results.bed`.")

# **Fix column misalignment by keeping everything after Frame in `Attribute` column**
expected_cols = 13  # Expected number of columns before "Attribute"
if intersected_df.shape[1] > expected_cols:
    intersected_df.iloc[:, expected_cols - 1] = intersected_df.iloc[:, expected_cols - 1:].apply(lambda x: ' '.join(x.dropna()), axis=1)
    intersected_df = intersected_df.iloc[:, :expected_cols]  # Trim extra columns

# **Assign corrected column names**
column_names = ["CHROM", "Start_pqtl", "End_pqtl", "ID", "GTF_CHROM", "Source", "Feature", 
                "GTF_Start", "GTF_End", "Score", "Strand", "Frame", "Attribute"]
intersected_df.columns = column_names

# Debugging: Print first rows to verify correct column mapping
print("🔹 Sample of corrected intersected_df:")
print(intersected_df.head())

# **Check if ID column is correctly extracted**
if "ID" not in intersected_df.columns or intersected_df["ID"].str.contains("HAVANA").all():
    raise ValueError(
        "❌ Error: 'ID' column is incorrect! It appears column misalignment has occurred in `intersected_results.bed`.\n"
        "👉 Run the following command in the terminal to check column alignment:\n"
        "   head /Users/shivank/Desktop/Hackathon/test/intersected_results.bed | cat -t"
    )

# Load Full pqtl data
pqtl_df = pd.read_csv(pqtl_combined_file, sep="\t", dtype=str)

# Debugging: Check column names before merging
print("🔹 pqtl_df Columns:", pqtl_df.columns.tolist())
print("🔹 intersected_df Columns:", intersected_df.columns.tolist())

# Ensure ID formatting is consistent before merging
pqtl_df["ID"] = pqtl_df["ID"].astype(str).str.strip()
intersected_df["ID"] = intersected_df["ID"].astype(str).str.strip()

# Debugging: Check if ID values match
print("🔹 Sample pqtl_df IDs:", pqtl_df["ID"].head().tolist())
print("🔹 Sample intersected_df IDs:", intersected_df["ID"].head().tolist())

# Merge pqtl with intersected data based on ID
merged_df = pqtl_df.merge(intersected_df, on="ID", how="inner")  # Keep all columns

# Debugging: Check if merge worked
print(f"✅ Merged dataset: {len(merged_df)} variants (before filtering)")

# Function to extract Gene_ID from the attribute column
def extract_gene_id(attribute_str):
    match = re.search(r'gene_id "([^"]+)"', attribute_str)  # Match `gene_id "ENSG00000184731.6"`
    return match.group(1).split(".")[0] if match else None  # Remove `.number`

# Extract Gene_ID **only after merging**
merged_df["Gene_ID"] = merged_df["Attribute"].apply(extract_gene_id)

# Load Gene List and remove `.number`
gene_list_df = pd.read_csv(gene_list_file, sep="\t", names=["Gene_ID"], dtype=str)
gene_list_df["Gene_ID"] = gene_list_df["Gene_ID"].astype(str).str.split(".").str[0]

# Convert gene list to a set for faster filtering
gene_set = set(gene_list_df["Gene_ID"])

# Filter merged data based on Gene_ID
filtered_df = merged_df[merged_df["Gene_ID"].isin(gene_set)]

print(f"✅ Final filtered pqtl data saved: {len(filtered_df)} variants (after filtering)")

# Save the final filtered dataset
filtered_df.to_csv(output_file, sep="\t", index=False)

print(f"✅ Final file saved: {output_file}")
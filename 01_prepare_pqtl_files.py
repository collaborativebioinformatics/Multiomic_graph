import pandas as pd
import glob
import re

### USER INPUTS ###
pqtl_files_path = "/test/pqtl_test_files/*.gz"
output_bed_file = "/test/pqtl_data.bed"
output_tsv_file = "/test/pqtl_combined.tsv"

# Get all pqtl files
pqtl_files = glob.glob(pqtl_files_path)

# Expected column structure
expected_columns = ["CHROM", "GENPOS", "ID", "ALLELE0", "ALLELE1", "A1FREQ", "INFO", "N", "TEST", "BETA", "SE", "CHISQ", "LOG10P", "EXTRA"]

# Read and concatenate pqtl files
pqtl_df_list = []
for file in pqtl_files:
    df = pd.read_csv(file, compression="gzip", sep="\s+", names=expected_columns, header=0, engine="python")

    # Drop duplicate headers appearing as rows
    if df.columns[0] in df.iloc[0].values:
        df = df.iloc[1:]

    pqtl_df_list.append(df)

# Concatenate all files
pqtl_df = pd.concat(pqtl_df_list, ignore_index=True)

# Extract POSITION from ID
def extract_position(id_str):
    match = re.match(r'[^:]+:(\d+):.*', str(id_str))
    return int(match.group(1)) if match else None

pqtl_df["Start"] = pqtl_df["ID"].apply(extract_position) - 1  # Convert to 0-based BED format
pqtl_df["End"] = pqtl_df["Start"] + 1  # BED format requires End = Start + 1

# Standardize chromosome format (Ensure consistency with GTF)
pqtl_df["CHROM"] = "chr" + pqtl_df["CHROM"].astype(str)

# Ensure columns are correctly formatted
pqtl_df["Start"] = pqtl_df["Start"].astype(int)
pqtl_df["End"] = pqtl_df["End"].astype(int)

# Sorting function that ensures CHROM is always a string
def chrom_sort_key(chrom):
    chrom = str(chrom)  # Convert to string (handles cases where CHROM is int)
    chrom = chrom.replace("chr", "")  # Remove "chr" prefix if present
    try:
        return int(chrom)  # Sort numerically (e.g., 1, 2, 3…)
    except ValueError:
        return {"X": 23, "Y": 24, "MT": 25}.get(chrom, float("inf"))  # Sort X, Y, MT last

# Ensure CHROM is a string before sorting
pqtl_df["CHROM"] = pqtl_df["CHROM"].astype(str)

# Sort by CHROM, Start, End
pqtl_df = pqtl_df.sort_values(by=["CHROM", "Start", "End"], key=lambda x: x.map(chrom_sort_key))

# Save BED file (for bedtools)
pqtl_bed = pqtl_df[["CHROM", "Start", "End", "ID"]]
pqtl_bed.to_csv(output_bed_file, sep="\t", index=False, header=False)

# Save full pqtl data file
pqtl_df.to_csv(output_tsv_file, sep="\t", index=False)

print(f"✅ pqtl_data.bed created with {len(pqtl_bed)} variants (for bedtools)")
print(f"✅ pqtl_combined.tsv created with full pqtl dataset: {len(pqtl_df)} rows")
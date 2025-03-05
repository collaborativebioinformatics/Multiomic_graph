import pandas as pd
import re

### USER INPUTS ###
gtf_file = "test/gencode.v47.annotation.gtf"
output_file = "test/gtf_data.bed"

# Read GTF file (keeping all columns)
gtf_columns = ["CHROM", "source", "feature", "start", "end", "score", "strand", "frame", "attribute"]
gtf_df = pd.read_csv(gtf_file, sep="\t", comment="#", names=gtf_columns)

# Extract gene_id and gene_type from the attribute column
def extract_gene_id(attribute_str):
    match = re.search(r'gene_id "([^"]+)"', attribute_str)
    return match.group(1) if match else None

def extract_gene_type(attribute_str):
    match = re.search(r'gene_type "([^"]+)"', attribute_str)
    return match.group(1) if match else None

gtf_df["Gene_ID"] = gtf_df["attribute"].apply(extract_gene_id)
gtf_df["Gene_Type"] = gtf_df["attribute"].apply(extract_gene_type)

# Filter for feature == "gene" and gene_type == "protein_coding" # & (gtf_df["Gene_Type"] == "protein_coding")
gtf_df = gtf_df[(gtf_df["feature"] == "gene") ]

# Convert to BED format (Start must be 0-based)
gtf_df["Start"] = gtf_df["start"] - 1  # Convert 1-based GTF format to 0-based BED

# Ensure `CHROM` format is consistent with pqtl_data.bed
gtf_df["CHROM"] = gtf_df["CHROM"].astype(str)
if not gtf_df["CHROM"].str.startswith("chr").all():
    gtf_df["CHROM"] = "chr" + gtf_df["CHROM"]

# Ensure integer formatting for BED format
gtf_df["Start"] = gtf_df["Start"].astype(int)
gtf_df["end"] = gtf_df["end"].astype(int)

# Save BED file with TAB separation (keeping all columns)
gtf_df.to_csv(output_file, sep="\t", index=False, header=False)

print(f"✅ gtf_data.bed created with {len(gtf_df)} genes")
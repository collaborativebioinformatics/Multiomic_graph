import subprocess
import pandas as pd

### USER INPUTS ###
pqtl_bed_file = "test/pqtl_data.bed"
gtf_bed_file = "test/gtf_data.bed"
output_file = "test/intersected_results.bed"

print("🔹 Running bedtools intersect...")

# Run BEDTools Intersect
bedtools_command = f"bedtools intersect -a {pqtl_bed_file} -b {gtf_bed_file} -wa -wb > {output_file}"
subprocess.run(bedtools_command, shell=True, check=True)

# Load the intersected results and remove empty rows if any
intersected_df = pd.read_csv(output_file, sep="\t", header=None)

# Remove any empty lines (if any exist)
intersected_df = intersected_df.dropna()

# Save cleaned intersected results
intersected_df.to_csv(output_file, sep="\t", index=False, header=False)

print(f"✅ Cleaned intersected results saved: {len(intersected_df)} rows remaining")
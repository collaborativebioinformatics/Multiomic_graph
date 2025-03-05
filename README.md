# Multiomic_graph

Group Members: 
Arth Banka,
 Siddharth Sabata,
 Shivank Sadasivan,
 Rachael Oluwakamiye abolade,
 Lars W. Ericson

Population-Specific Multi-omics Graph Generation for a Target Protein
This repository contains the code and resources for generating population-specific multi-omics gene graphs that integrate genetic variants and protein expression data. The project aims to link genomic variation with protein abundance using graph-based representations, facilitating advanced analytics and precision medicine applications.

Overview
Genetic variants can influence protein expression, leading to population-specific differences in biological processes and disease susceptibility. This project integrates protein quantitative trait loci (pQTL) data with population-specific genetic variants to construct gene-specific graphs. These graphs represent both reference genomic sequences and variant information, enabling downstream analysis through graph-based machine learning models.

Why Graphs?
Graphs enable advanced graph-based analytics and machine learning applications.

They facilitate precision medicine by linking genomic variation with protein expression.

Graphs provide a foundation for developing predictive models of variant impacts on protein expression.

Key Features
Variant Nodes:

Global genomic position

Reference and alternate allele sequences

Effect size (beta) of the variant

-log10 of the p-value

Direction of effect

Reference Nodes:

Global genomic position

Local index

Nucleotide sequence



Methods:
1. Target Protein pQTL Data Processing
The workflow begins with the processing of protein quantitative trait loci (pQTL) data for the target protein. This dataset consists of chromosome-specific variant files, one for each of the 23 chromosomes. These files are used to identify genetic variants associated with the target protein. The variants are extracted and organized by chromosome to facilitate downstream filtering and analysis. This step ensures that all relevant genetic information for the target protein is included in the pipeline.

2. Filtering Variants Based on User Input
The extracted pQTL data is filtered based on user-defined criteria, which currently involves parsing the INFO column of the pQTL dataset. The filtering criteria can be customized to include specific variant properties, such as allele frequency, effect size, or functional annotations. This step reduces the dataset to variants that meet user-specified thresholds, ensuring that only biologically or clinically relevant variants are included in subsequent analyses.

3. Gene Annotations and Variant Mapping
Gene annotation data in GTF (General Transfer Format) is used to map filtered variants to their corresponding genomic features. For each variant, its position is checked against annotated regions such as genes, untranslated regions (UTRs), and exons. If a variant falls within one of these regions, additional information such as the Ensembl gene ID is retrieved and linked to the variant. The current implementation focuses on protein-coding genes, but this can be extended to other gene types based on user input. This step ensures accurate mapping of variants to functional genomic elements.

4. Gene Sequence Extraction
The reference sequence for each gene identified in the previous step is retrieved using its Ensembl gene ID. This sequence forms the basis for constructing a linear "gene graph." The gene graph represents the reference sequence as a series of connected nodes, where each node corresponds to a segment of the sequence. This step provides a foundational structure for integrating variant information into the graph.

5. Variant Integration into Gene Graphs
Variants mapped to genes are integrated into their respective gene graphs by modifying or adding nodes and edges. Each node in the graph stores metadata such as nucleotide sequence, length, position, strand orientation, chromosome number, and genomic feature type (e.g., gene, exon, UTR). For variant nodes, additional properties such as base effect (e.g., SNPs, insertions, deletions) and statistical significance (e.g., logP value) are also stored. Reference sequence nodes do not include these additional properties. All variants are connected to both upstream and downstream nodes in the graph to ensure continuity.

6. Construction of Target Protein Genome Graph
The final step involves combining all gene graphs associated with the target protein into a comprehensive genome graph. This graph incorporates both reference sequences and variant information for all relevant genes across chromosomes. The resulting structure provides a holistic view of genetic variation affecting the target protein and enables downstream analyses such as path traversal or visualization of alternative haplotypes.

7. Parallelization
To optimize performance, especially when processing large datasets or multiple chromosomes simultaneously, parallelization is explored within this workflow. Computationally intensive steps such as filtering variants or constructing gene graphs can be parallelized across multiple processors or distributed computing environments. This ensures scalability and efficiency when handling high-throughput sequencing data.



![image](https://github.com/user-attachments/assets/3a81935f-a8e2-434f-8589-bae07b105f82)





Here’s the updated **GitHub README** based on your workflow, which consists of four scripts. The workflow section, installation, next steps, and results have been revised accordingly.

---

# Population-Specific Multi-omics Graph Generation for a Target Protein

This repository contains the code and resources for generating **population-specific multi-omics gene graphs** that integrate genetic variants and protein expression data. The project links genomic variation with protein abundance using graph-based representations, facilitating advanced analytics and precision medicine applications.

---

## Overview

Genetic variants can influence protein expression, leading to population-specific differences in biological processes and disease susceptibility. This project integrates **protein quantitative trait loci (pQTL)** data with population-specific genetic variants to construct **gene-specific graphs**. These graphs represent both reference genomic sequences and variant information, enabling downstream analysis through graph-based machine learning models.

---

## Workflow

The pipeline consists of **four scripts**, each performing a specific task in the workflow:

### Script 1: pQTL Data Processing
This script processes multiple pQTL files, concatenates them into a single dataset, extracts genomic positions from variant IDs, and converts the data into BED format for compatibility with downstream tools like BEDTools. The output includes:
- A combined pQTL dataset (`pqtl_combined.tsv`) with full variant information.
- A BED file (`pqtl_data.bed`) containing chromosome, start, end positions, and variant IDs for intersection analysis.

### Script 2: GTF File Processing
This script processes a GTF annotation file to extract gene regions and convert them into BED format. It filters for protein-coding genes and ensures consistency in chromosome naming conventions. The output is a BED file (`gtf_data.bed`) containing gene coordinates for intersection analysis.

### Script 3: Variant-Gene Intersection
This script uses `bedtools intersect` to compute overlaps between pQTL variants and gene regions. It identifies variants that fall within annotated genes and outputs intersected results in BED format (`intersected_results.bed`).

### Script 4: Final Filtering
This script merges the intersected results with the full pQTL dataset, extracts gene IDs from the annotation column, and filters variants based on a user-provided list of target genes (`test_gene_list.csv`). The output is a filtered dataset (`final_filtered_pqtl.tsv`) containing only variants associated with target genes.

---

## Installation

### Prerequisites

- Python >= 3.8
- Required Python libraries: `pandas`, `numpy`, `re`
- BEDTools (command-line tool)

### Setup Instructions

1. Clone this repository:
   ```bash
   git clone https://github.com/your-repo-name/genomic-graph-generation.git
   cd genomic-graph-generation
   ```

2. Install Python dependencies:
   ```bash
   pip install pandas numpy
   ```

3. Install BEDTools:
   Follow instructions from [BEDTools documentation](https://bedtools.readthedocs.io/en/latest/content/installation.html).

4. Ensure access to required input files:
   - pQTL dataset files (`*.gz`)
   - GTF annotation file (e.g., `gencode.v47.annotation.gtf`)
   - Target gene list (e.g., `test_gene_list.csv`)

---

## Usage

Run each script sequentially to complete the pipeline:

### Step 1: Process pQTL Files
```bash
python script1_pqtl_processing.py
```
Input:
- Path to pQTL files (`*.gz`)
Output:
- Combined pQTL dataset (`pqtl_combined.tsv`)
- BED file for variants (`pqtl_data.bed`)

---

### Step 2: Process GTF File
```bash
python script2_gtf_processing.py
```
Input:
- Path to GTF file (e.g., `gencode.v47.annotation.gtf`)
Output:
- BED file for gene regions (`gtf_data.bed`)

---

### Step 3: Compute Variant-Gene Intersection
```bash
python script3_bedtools_intersect.py
```
Input:
- Variant BED file (`pqtl_data.bed`)
- Gene BED file (`gtf_data.bed`)
Output:
- Intersected results (`intersected_results.bed`)

---

### Step 4: Filter Variants by Target Genes
```bash
python script4_final_filtering.py
```
Input:
- Intersected results (`intersected_results.bed`)
- Combined pQTL dataset (`pqtl_combined.tsv`)
- Target gene list (`test_gene_list.csv`)
Output:
- Final filtered dataset (`final_filtered_pqtl.tsv`)

---

## Results

The pipeline generates the following outputs:

1. **Intermediate Outputs**:
   - `pqtl_combined.tsv`: Full pQTL dataset.
   - `pqtl_data.bed`: BED file for variant positions.
   - `gtf_data.bed`: BED file for gene regions.
   - `intersected_results.bed`: Variants intersecting gene regions.

2. **Final Output**:
   - `final_filtered_pqtl.tsv`: Filtered dataset containing only variants associated with target genes.

Example log from Script 4:
```
✅ Merged dataset: 3,237 variants (before filtering)
✅ Final filtered pqtl data saved: 1,245 variants (after filtering)
✅ Final file saved: final_filtered_pqtl.tsv
```

---

## Next Steps

1. **Graph Construction**:
   - Use the filtered pQTL data to construct population-specific gene graphs.
   - Represent reference sequences and variant nodes using graph-based libraries like [PyTorch Geometric](https://pytorch-geometric.readthedocs.io/).

2. **Machine Learning Integration**:
   - Train graph neural networks (GNNs) to predict variant impacts on protein expression.
   - Explore clustering or classification tasks using graph embeddings.

3. **Visualization**:
   - Develop interactive tools for visualizing gene graphs and exploring variant effects.

4. **Scalability**:
   - Optimize pipeline performance for large datasets.
   - Parallelize computationally intensive steps (e.g., intersection analysis).

---

## Contributors

This project was developed as part of a hackathon by:

- Siddharth Sabata  
- Shivank Sadasivan  
- Lars Warren Ericson  
- Rachael Oluwakamiye Abolade  
- Arth Banka  

---

## License

This project is licensed under the MIT License. See the [LICENSE](LICENSE) file for details.

---

## Acknowledgments

We thank the organizers of the hackathon for providing a collaborative environment and resources to work on this project.

---

This README now reflects your workflow accurately while maintaining clarity and professionalism. Let me know if you'd like further refinements!

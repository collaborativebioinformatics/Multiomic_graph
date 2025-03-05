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



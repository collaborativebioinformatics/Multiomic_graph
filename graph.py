#!/usr/bin/env python3
"""
Genome Graph Generator

This script builds PyTorch Geometric graph representations of genes with their variants.
Each graph represents a gene's reference sequence as a backbone with variant nodes.

Usage:
    python graph.py --csv variants.csv --email your.email@example.com

Required CSV format:
    The input CSV must contain the following columns:
    - GENEID: Ensembl gene ID (e.g., ENSG00000092847.14)
    - ALLELE0: Reference allele sequence
    - ALLELE1: Alternate allele sequence
    - POSITION: Genomic position (chromosome coordinate)
    - BETA: Effect size
    - LOG10P: -log10(p-value)
    
    Optional columns:
    - direction: Strand direction ("+" or "-"), defaults to "+"

Outputs:
    - A text file with the list of gene IDs in order (saved in outputs directory)
    - A compressed PyTorch file (.pt.gz) containing the list of PyG Data objects (saved in outputs directory)
    
Requirements:
    - pandas
    - requests
    - biopython
    - torch
    - torch_geometric
"""

import pandas as pd
import requests
from Bio import Entrez
import torch
from torch_geometric.data import Data
import math
import argparse
import sys
import os
import gzip

############################
# PART A: Fetching Reference from Ensembl
############################

def get_genomic_coordinates(ensembl_id):
    """Retrieve chromosome, start, and end positions from Ensembl REST API."""
    base_id = ensembl_id.split('.')[0]  # Remove version number if present
    url = f"https://rest.ensembl.org/lookup/id/{base_id}?expand=0"
    
    response = requests.get(url, headers={"Content-Type": "application/json"})
    response.raise_for_status()
    
    data = response.json()
    return {
        'chrom': data['seq_region_name'],
        'start': data['start'] - 1,  # Convert to 0-based
        'end': data['end']
    }

def get_ncbi_refseq_id(chrom):
    """Construct NCBI RefSeq ID from chromosome number (for GRCh38)."""
    # For a real pipeline, you may need a more robust mapping from Ensembl chromosome to NCBI.
    # Below is a simple placeholder that handles numeric, 'X', 'Y'.
    if chrom in ['X', 'Y']:
        # For GRCh38, X -> NC_000023.11, Y -> NC_000024.10
        # Hard-code or create a small dictionary if you prefer:
        return "NC_000023.11" if chrom == 'X' else "NC_000024.10"
    else:
        # For autosomes, e.g. '1' -> "NC_000001.11", '2' -> "NC_000002.12", etc.
        # You might need a dictionary mapping for each chromosome version.
        # For demonstration, we'll just do "NC_00000{chrom}.11" if chrom is 1..9
        # This is not guaranteed correct for all chromosomes, but shows the idea.
        num = int(chrom)  # if chrom is e.g. '3'
        return f"NC_{num:06d}.11"  # e.g. "NC_000003.11"

def fetch_fasta(ensembl_id, email="youremail@example.com"):
    """Fetch FASTA sequence from NCBI using Ensembl ID coordinates."""
    coordinates = get_genomic_coordinates(ensembl_id)
    chrom = coordinates['chrom']
    start = coordinates['start']
    end = coordinates['end']
    
    refseq_id = get_ncbi_refseq_id(chrom)
    
    Entrez.email = email  # NCBI requirement
    handle = Entrez.efetch(
        db="nuccore",
        id=refseq_id,
        rettype="fasta",
        retmode="text",
        seq_start=start + 1,  # efetch expects 1-based
        seq_stop=end
    )
    fasta_str = handle.read()
    handle.close()
    
    # The FASTA string includes header lines plus sequence lines.
    # We'll parse out just the sequence as a continuous string:
    lines = fasta_str.split('\n')
    seq_lines = [ln for ln in lines if not ln.startswith('>')]
    reference_seq = "".join(seq_lines).upper()
    
    return reference_seq, start  # returning the 0-based start

############################
# PART B: PyG Graph Construction
############################

def create_reference_graph_pyg(reference, gene, start_index=0):
    """
    Build a reference backbone for the genome as an intermediate dictionary.
    Each base in the reference becomes a node with attributes:
      - position: global genomic coordinate = start_index + i
      - local_index: i
      - nucleotide: the base (str)
      - type: "reference"
      - gene: gene name
    """
    nodes = {}
    edges = []
    n = len(reference)
    
    for i, base in enumerate(reference):
        global_pos = start_index + i
        node_id = str(i)
        nodes[node_id] = {
            "position": global_pos,
            "local_index": i,
            "nucleotide": base,
            "type": "reference",
            "gene": gene
        }
    for i in range(n - 1):
        edges.append((str(i), str(i+1)))
    
    return nodes, edges

def add_variant(nodes, edges, reference,
                global_p, ref_seq, alt_seq,
                gene, beta=None, log10_p=None,
                direction="+",
                start_index=0):
    """
    Add one variant bubble node for a mutation at global position `global_p`.
    local_p = global_p - start_index
    """
    n = len(reference)
    local_p = global_p - start_index
    
    if local_p < 0 or local_p >= n:
        print(f"Warning: variant position {global_p} (local={local_p}) out of range.")
        return
    
    var_node_id = f"VAR_{global_p}_{ref_seq}->{alt_seq}"
    nodes[var_node_id] = {
        "position": global_p,
        "ref": ref_seq,
        "alt": alt_seq,
        "alt_length": len(alt_seq),
        "ref_length": len(ref_seq),
        "beta": beta,
        "log10_p": log10_p,
        "direction": direction,
        "type": "variant",
        "gene": gene
    }
    
    # left anchor
    if local_p > 0:
        edges.append((str(local_p-1), var_node_id))
    
    # right anchor
    end_anchor = local_p + len(ref_seq)
    if end_anchor < n:
        edges.append((var_node_id, str(end_anchor)))

def add_variants_to_graph(nodes, edges,
                          reference, variants, gene,
                          start_index=0):
    """
    For each variant in `variants`, call add_variant.
    Each variant dict has:
      - position (global)
      - ref, alt
      - beta, log10_p
      - direction
    """
    for var in variants:
        add_variant(
            nodes, edges, reference,
            global_p=var['position'],
            ref_seq=var['ref'],
            alt_seq=var['alt'],
            gene=gene,
            beta=var.get("beta", None),
            log10_p=var.get("log10_p", None),
            direction=var.get("direction", "+"),
            start_index=start_index
        )

def create_genome_graph_pyg(reference, variants, gene, start_index=0):
    """
    1) Create reference backbone (undirected).
    2) Add variant bubble nodes.
    3) Convert to PyG Data object.
    """
    nodes, edges = create_reference_graph_pyg(reference, gene, start_index=start_index)
    add_variants_to_graph(nodes, edges, reference, variants, gene, start_index=start_index)
    
    node_ids = list(nodes.keys())
    node2index = {nid: idx for idx, nid in enumerate(node_ids)}
    
    features = []
    node_labels = {}
    for nid in node_ids:
        attr = nodes[nid]
        features.append([attr["position"]])
        if attr.get("type") == "reference":
            node_labels[node2index[nid]] = attr.get("nucleotide", "")
        else:
            node_labels[node2index[nid]] = attr.get("alt", "")
    
    x = torch.tensor(features, dtype=torch.float)
    
    edge_index_list = []
    for src, tgt in edges:
        src_idx = node2index[src]
        tgt_idx = node2index[tgt]
        edge_index_list.append([src_idx, tgt_idx])
        edge_index_list.append([tgt_idx, src_idx])
    edge_index = torch.tensor(edge_index_list, dtype=torch.long).t().contiguous()
    
    data = Data(x=x, edge_index=edge_index)
    data.node_info = nodes
    data.node_labels = node_labels
    data.index_to_nid = {node2index[nid]: nid for nid in node_ids}
    
    return data

############################
# PART C: Convert DF Rows -> Variants
############################

def df_to_variants(df):
    """
    Convert each row in df into a dict with:
      - position
      - ref
      - alt
      - beta
      - log10_p
      - direction
      - gene (optional, but we'll store it in the big script below)
    """
    # We'll let the grouping function handle the gene logic.
    # This function only extracts the needed columns.
    variant_dicts = []
    for _, row in df.iterrows():
        variant = {
            "position": row["POSITION"],
            "ref": row["ALLELE0"],
            "alt": row["ALLELE1"],
            "beta": row["BETA"],
            "log10_p": row["LOG10P"],
            "direction": row.get("direction", "+"),
        }
        variant_dicts.append(variant)
    return variant_dicts

############################
# PART D: Master Function to Build Graphs for Each Gene
############################

def build_genome_graphs_for_all_genes(df, output_textfile="genes_order.txt", email="youremail@example.com"):
    """
    1. Group the DataFrame by GENEID (which is an Ensembl ID).
    2. For each gene:
       - Fetch the reference sequence from Ensembl/NCBI (fasta).
       - Parse the variants for that gene (df_to_variants).
       - Create a PyG graph with create_genome_graph_pyg.
       - Append the graph to a list.
    3. Write the Ensembl IDs in order to output_textfile.
    4. Return the list of PyG graphs.
    """
    # We'll store final results here
    graph_list = []
    gene_list = []  # keep track of gene IDs in order
    
    # Group by 'GENEID'
    grouped = df.groupby("GENEID")
    total_genes = len(grouped)
    
    print(f"Processing {total_genes} genes...")
    
    for i, (gene_id, sub_df) in enumerate(grouped):
        progress = (i + 1) / total_genes * 100
        print(f"[{progress:.1f}%] Processing gene {i+1}/{total_genes}: {gene_id}")
        
        # 1) fetch the reference
        #    - We'll assume gene_id is an Ensembl ID (e.g., "ENSG00000092847.14")
        #    - fetch_fasta returns (reference_seq, start_coord)
        try:
            reference_seq, start_coord = fetch_fasta(gene_id, email=email)
        except Exception as e:
            print(f"  Error fetching FASTA for gene {gene_id}: {e}")
            continue
            
        print(f"  Retrieved reference sequence ({len(reference_seq)} bp) starting at position {start_coord}")
        
        # 2) parse the variants for that gene
        #    - We'll convert the sub_df to a list of variant dicts
        variants = []
        for _, row in sub_df.iterrows():
            variant = {
                "position": row["POSITION"],  # global coordinate
                "ref": row["ALLELE0"],
                "alt": row["ALLELE1"],
                "beta": row["BETA"],
                "log10_p": row["LOG10P"],
                "direction": row.get("direction", "+"),
            }
            variants.append(variant)
            
        print(f"  Found {len(variants)} variants for this gene")
        
        # 3) create the PyG graph
        data = create_genome_graph_pyg(
            reference_seq,
            variants,
            gene=gene_id,
            start_index=start_coord  # so that local_p = position - start_coord
        )
        
        # 4) store in our lists
        graph_list.append(data)
        gene_list.append(gene_id)
    
    # 5) Write genes to file
    with open(output_textfile, "w") as f:
        for g in gene_list:
            f.write(f"{g}\n")
    
    return graph_list

def save_genome_graphs(graphs, output_file):
    """
    Save a list of PyTorch Geometric Data objects to a compressed file.
    
    Args:
        graphs: List of PyTorch Geometric Data objects
        output_file: Path to save the compressed graphs
    """
    print(f"Saving {len(graphs)} graphs to {output_file}...")
    try:
        with gzip.open(output_file, "wb") as f:
            torch.save(graphs, f)
        return True
    except Exception as e:
        print(f"Error saving graphs: {e}")
        return False

############################
# Example Usage
############################
if __name__ == "__main__":
    # Set up command-line argument parsing
    parser = argparse.ArgumentParser(
        description="Build genome graphs from a CSV file of variant data",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    
    parser.add_argument(
        "--csv", "-c", 
        required=True,
        help="Path to the CSV file containing variant data"
    )
    
    parser.add_argument(
        "--email", "-e", 
        required=True,
        help="Email address for NCBI services"
    )
    
    parser.add_argument(
        "--output-dir", "-d",
        default="outputs",
        help="Directory to store output files"
    )
    
    parser.add_argument(
        "--gene-list", "-l",
        default="genes_order.txt",
        help="Name of the file to store gene order (will be saved in output directory)"
    )
    
    parser.add_argument(
        "--graphs-output", "-g",
        default="genome_graphs.pt.gz",
        help="Name of the file to store the compressed PyTorch Geometric graphs (will be saved in output directory)"
    )
    
    parser.add_argument(
        "--verbose", "-v", 
        action="store_true",
        help="Print detailed progress information"
    )
    
    # Parse arguments
    args = parser.parse_args()
    
    # Validate CSV file exists
    if not os.path.isfile(args.csv):
        print(f"Error: CSV file '{args.csv}' not found.")
        sys.exit(1)
    
    # Create output directory if it doesn't exist
    if not os.path.exists(args.output_dir):
        print(f"Creating output directory: {args.output_dir}")
        os.makedirs(args.output_dir)
    
    # Prepare output file paths
    gene_list_path = os.path.join(args.output_dir, args.gene_list)
    graphs_output_path = os.path.join(args.output_dir, args.graphs_output)
    
    # Load CSV into DataFrame
    try:
        print(f"Loading variant data from {args.csv}...")
        df = pd.read_csv(args.csv)
        
        # Check required columns
        required_columns = ["GENEID", "ALLELE0", "ALLELE1", "POSITION", "BETA", "LOG10P"]
        missing_columns = [col for col in required_columns if col not in df.columns]
        
        if missing_columns:
            print(f"Error: Missing required columns in CSV: {', '.join(missing_columns)}")
            print(f"Required columns are: {', '.join(required_columns)}")
            sys.exit(1)
            
        print(f"Loaded {len(df)} variants across {df['GENEID'].nunique()} genes.")
        
    except Exception as e:
        print(f"Error loading CSV file: {e}")
        sys.exit(1)
    
    # Run the main function
    print(f"Building genome graphs...")
    try:
        all_graphs = build_genome_graphs_for_all_genes(
            df, 
            output_textfile=gene_list_path, 
            email=args.email
        )
        
        print(f"Successfully built {len(all_graphs)} graphs.")
        for i, gdata in enumerate(all_graphs):
            if args.verbose:
                print(f"Graph {i}: gene={list(gdata.node_info.values())[0].get('gene','?')}, nodes={gdata.num_nodes}, edges={gdata.num_edges}")
        
        print(f"Gene order saved to {gene_list_path}")
        
        # Save graphs to file
        if save_genome_graphs(all_graphs, graphs_output_path):
            print(f"Graphs saved to {graphs_output_path}")
    
    except Exception as e:
        print(f"Error building genome graphs: {e}")
        sys.exit(1)

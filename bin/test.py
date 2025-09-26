#!/usr/bin/env python
"""
RefSeq Genome Size Extractor

This script downloads and parses RefSeq assembly summary files from NCBI
for specified taxonomic groups, calculates genome sizes for each assembly,
and outputs both per-assembly and per-species average genome size tables.

Author: Xiaoli Dong
Version: 1.0
Date: 2025-07-18

Usage example:
  python genome_size_estimator.py -g bacteria -o1 assemblies.tsv -o2 organisms.tsv
"""

import pandas as pd
import requests
from io import StringIO
import argparse

def download_and_parse_assembly_summary(url):
    response = requests.get(url)
    response.raise_for_status()
    header_line = next(line for line in response.text.splitlines()
                       if line.startswith("#assembly_accession"))
    columns = header_line.lstrip("#").split("\t")
    lines = [line for line in response.text.splitlines()
             if not line.startswith("#")]
    data_str = "\n".join(lines)
    df = pd.read_csv(StringIO(data_str), sep="\t", header=None, names=columns)
    return df

def main():
    parser = argparse.ArgumentParser(
        description="Download and parse NCBI RefSeq assembly summaries.")
    parser.add_argument(
        "-g", "--groups",
        nargs="+",
        default=["all"],
        help="Taxonomic groups to download. Options: archaea bacteria fungi viral plant protozoa invertebrate vertebrate_mammalian vertebrate_other all"
    )
    parser.add_argument(
        "-o1", "--output_assemblies",
        default="assemblies_with_genome_size.tsv",
        help="Output file for individual assemblies (default: assemblies_with_genome_size.tsv)"
    )
    parser.add_argument(
        "-o2", "--output_organisms",
        default="organism_average_genome_size.tsv",
        help="Output file for average genome size per organism (default: organism_average_genome_size.tsv)"
    )
    args = parser.parse_args()

    all_groups = [
        "archaea", "bacteria", "fungi", "viral", "plant",
        "protozoa", "invertebrate", "vertebrate_mammalian", "vertebrate_other"
    ]

    if "all" in args.groups:
        groups = all_groups
    else:
        groups = [g.lower() for g in args.groups]
        invalid = [g for g in groups if g not in all_groups]
        if invalid:
            print(f"Error: invalid group(s) specified: {invalid}")
            return

    base_url = "https://ftp.ncbi.nlm.nih.gov/genomes/refseq"
    df_list = []

    for group in groups:
        url = f"{base_url}/{group}/assembly_summary.txt"
        print(f"Downloading {group} summary...")
        try:
            df_group = download_and_parse_assembly_summary(url)
            df_group['group'] = group
            df_list.append(df_group)
        except Exception as e:
            print(f"Failed to download or parse {group}: {e}")

    if not df_list:
        print("No data downloaded.")
        return

    df_all = pd.concat(df_list, ignore_index=True)
    df_all['genome_size'] = pd.to_numeric(
        df_all['genome_size'], errors='coerce')
    df_latest = df_all[df_all["version_status"] == "latest"]
    df_valid = df_latest.dropna(subset=['genome_size'])

    # Prepare assembly output table
    output_df = df_valid[[
        'assembly_accession', 'taxid', 'species_taxid', 'organism_name',
        'infraspecific_name', 'assembly_level', 'version_status', 'genome_size'
    ]]
    output_df = output_df.sort_values(
        by=['organism_name', 'genome_size'], ascending=[True, False])

    output_df.to_csv(args.output_assemblies, index=False, sep=',')
    print(f"Saved individual assembly records to: {args.output_assemblies}")

    # Separate complete genomes and others
    complete_df = df_valid[df_valid["assembly_level"] == "Complete Genome"]
    others_df = df_valid[df_valid["assembly_level"] != "Complete Genome"]

    # Find species_taxids with at least one complete genome
    species_with_complete = complete_df["species_taxid"].unique()

    # Group complete genomes
    grouped_complete = (
        complete_df.groupby("species_taxid", as_index=False)
        .agg(
            organism_name=("organism_name", "first"),
            avg_genome_size=("genome_size", "mean"),
            n_assemblies=("genome_size", "count"),
            used_complete_only=("assembly_level", lambda x: True)
        )
    )

    # Group other genomes only for species not in complete set
    remaining_df = others_df[~others_df["species_taxid"].isin(species_with_complete)]
    grouped_other = (
        remaining_df.groupby("species_taxid", as_index=False)
        .agg(
            organism_name=("organism_name", "first"),
            avg_genome_size=("genome_size", "mean"),
            n_assemblies=("genome_size", "count"),
            used_complete_only=("assembly_level", lambda x: False)
        )
    )

    # Combine and round
    grouped_df = pd.concat([grouped_complete, grouped_other], ignore_index=True)
    grouped_df["avg_genome_size"] = grouped_df["avg_genome_size"].round(0).astype(int)

    # Sort by organism_name then genome size
    grouped_df = grouped_df.sort_values(
        by=['organism_name', 'avg_genome_size'], ascending=[True, False])

    # Save grouped summary
    grouped_df.to_csv(args.output_organisms, index=False, sep=',')
    print(f"Saved organism-level average genome sizes to: {args.output_organisms}")

    # Overall summary stats
    print(f"\nTotal entries with genome size: {len(df_valid)}")
    print(f"Average genome size (bp): {df_valid['genome_size'].mean():,.0f}")


if __name__ == "__main__":
    main()

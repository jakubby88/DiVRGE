#this code takes a timo or group of timo outputs and calculates dN/dS

import argparse
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Data import CodonTable
import pandas as pd

#genetic code
genetic_code = CodonTable.unambiguous_dna_by_id[1]

# function to count expected numbers of synonymous and nonsynonymous mutations at a base position
def calc_syn_nonsyn(codon, position):
    syn = 0
    nonsyn = 0
    
    bases = ['A', 'T', 'C', 'G']
    for base in bases:
        if base != codon[position - 1]:  # only consider mutations
            mutated_codon = list(codon)
            mutated_codon[position - 1] = base
            mutated_codon = ''.join(mutated_codon)
            
            if genetic_code.forward_table.get(mutated_codon, '*') == genetic_code.forward_table.get(codon, '*'):
                syn += 1
            else:
                nonsyn += 1
                
    # Return ratio of synonymous and nonsynonymous at a site bsed on 3 possible mutations
    return syn / 3, nonsyn / 3


def parse_genbank(genbank_file):
    # parse genbank file
    coding_regions = []
    with open(genbank_file, "r") as gb_file:
        for record in SeqIO.parse(gb_file, "genbank"):
            for feature in record.features:
                if feature.type == "CDS":
                    gene = feature.qualifiers.get("gene", ["unknown"])[0]
                    start = int(feature.location.start)
                    end = int(feature.location.end)
                    strand = feature.location.strand
                    coding_regions.append({"gene": gene, "start": start, "end": end, "strand": strand})
    return coding_regions
def calculate_sites(fasta_file, coding_regions):
    ref_genome = SeqIO.read(fasta_file, "fasta")
    results = []

    # process only the coding region genes
    for region in coding_regions:
        gene = region["gene"]
        start = region["start"]
        end = region["end"]
        
        for pos in range(start, end):
            relative_pos = pos - start
            
            # codon position and context
            if relative_pos % 3 == 0:
                codon_pos = 1
                codon_start = pos
                codon_end = pos + 3
            elif relative_pos % 3 == 2:
                codon_pos = 3
                codon_start = pos - 2
                codon_end = pos + 1
            elif relative_pos % 3 == 1:
                codon_pos = 2
                codon_start = pos - 1
                codon_end = pos + 2
            
            # codon from reference genome
            ref_codon = ref_genome.seq[codon_start:codon_end].upper()
            
            # skip if codon is incomplete
            if len(ref_codon) != 3:
                continue
            
            # calc synonymous and nonsynonymous counts
            syn_count, nonsyn_count = calc_syn_nonsyn(ref_codon, codon_pos)
            
            results.append({
                "gene": gene,
                "ntpos": pos + 1,
                "Synonymous_Sites": syn_count,
                "Nonsynonymous_Sites": nonsyn_count,
                "codon": ref_codon,
                "codon_position": codon_pos
            })

    results_df = pd.DataFrame(results)
    print("Expected sites calculated")
    return results_df


#take timo output and calculate dN/dS
def calculate_dn_ds_from_timo(input_csv, coding_regions, expected_values, additional_regions_csv=None):
    annotated_variants = pd.read_csv(input_csv)
    # If an additional regions CSV is provided, add these regions
    if additional_regions_csv:
        additional_regions = pd.read_csv(additional_regions_csv)
        for _, row in additional_regions.iterrows():
            coding_regions.append({
                "gene": row["gene"],
                "start": row["start"],
                "end": row["end"]
            })

    # Calculate dN/dS for each region
    dn_ds_results = []
    for region in coding_regions:
        gene, start, end = region["gene"], region["start"], region["end"]
        
        # Get expected and actual sites
        expected = expected_values.query(f"{start} <= ntpos <= {end}")
        actual = annotated_variants.query(f"{start} <= ntpos <= {end}")
        
        # Sum synonymous and nonsynonymous sites
        expected_syn = expected["Synonymous_Sites"].sum()
        expected_non = expected["Nonsynonymous_Sites"].sum()
        actual_syn = actual["synonymous"].sum()
        actual_non = actual["nonsynonymous"].sum()
        
        # Calculate dn/ds ratio
        if expected_syn > 0 and actual_syn > 0:
            dn_ds = (actual_non / actual_syn) / (expected_non / expected_syn)
        else:
            dn_ds = None
        
        dn_ds_results.append({"Gene": gene, "DN_DS": dn_ds})

    print("DN/DS calculation complete")
    return pd.DataFrame(dn_ds_results)

  

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Calculate number of nonsynonymous and synonymous sites in coding regions.")
    parser.add_argument("--output_file", type=str, required=True, help="Name of output CSV file")
    parser.add_argument("--ref_file", type=str, required=True, help="Name of reference FASTA file")
    parser.add_argument("--genbank_file", type=str, required=True, help="Path to the GenBank file")
    parser.add_argument("--input_file", type=str, required=True, help="Path to output timo file")
    parser.add_argument("--regions", type=str, help="Path to optional input csv that has regions of interest of DNDS. columns region, start (ntpos start), end(ntpos end). If no file 'None'")
    args = parser.parse_args()

    coding_regions = parse_genbank(args.genbank_file)
    expected_values = calculate_sites(args.ref_file, coding_regions)
    if args.regions is not None:
        dn_ds_results = calculate_dn_ds_from_timo(args.input_file, coding_regions, expected_values, args.regions)
    else:
        dn_ds_results = calculate_dn_ds_from_timo(args.input_file, coding_regions, expected_values)
    dn_ds_results.to_csv(args.output_file, index=False, header=True)
    print(f"dN/dS results saved to {args.output_file}")
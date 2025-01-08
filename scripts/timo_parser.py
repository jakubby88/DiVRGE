#this code takes a timo or group of timo outputs and calculates dN/dS
#certainly needs work to be more efficent, example more efficent use of ref fasta dict
#more importantly bc this is a "data anlysis" step variants being real needs to be a factor so need to build a p-value into timo

import argparse
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Data import CodonTable
import pandas as pd

#genetic code
genetic_code = CodonTable.unambiguous_dna_by_id[1]

def translate_codon(codon):
    codon = codon.upper()
    return genetic_code.forward_table.get(codon, "X")
def parse_fasta_to_dict(fasta_file):
    ref_fasta = {}
    for record in SeqIO.parse(fasta_file, "fasta"):
        for i, base in enumerate(record.seq):
            ref_fasta[i + 1] = base
    
    return ref_fasta

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

def parse_timo_output(timo_output, ref_fasta, major_threshold=0.5, minor_threshold=0.2):
    timo_data= pd.read_csv(timo_output)
    #for speed make the fasta file a dictionary 
    ref_fasta_dict = parse_fasta_to_dict(ref_fasta)

    filtered_variants = []
    for _, row in timo_data.iterrows():
        ntpos = row['ntpos']
        major_base = row['major']
        major_freq = row['majorfreq']
        minor_base = row['minor']
        minor_freq = row['minorfreq']
        
        # Get the reference base for this position
        ref_base = ref_fasta_dict.get(int(ntpos)) 
        if ref_base:
            # keep as concenous mutations, regardless of frequency
            if major_base != ref_base and major_base != minor_base and major_base != "N" and major_freq >= major_threshold:
                filtered_variants.append(row)
            
            # Keep minor bases only if they differ from the reference and have frequency >= 20%
            if minor_base != ref_base and minor_base != major_base and minor_base != "N" and minor_freq >= minor_threshold:
                filtered_variants.append(row)
    print(f"Number of filtered variants: {len(filtered_variants)}")
    return pd.DataFrame(filtered_variants)
def annotate_timo_output(filtered_variants, ref_genome, coding_regions):
    ref_genome = SeqIO.read(ref_genome, "fasta")
    nonsynonymous_list = []
    synonymous_list = []

    for _, row in filtered_variants.iterrows():
        ntpos = row['ntpos']
        major_base = row['major']

        # find the coding region containing the mutation
        gene_info = next((region for region in coding_regions if region['start'] <= ntpos <= region['end']), None)

        if gene_info:
            start = gene_info['start']

            relative_pos = ntpos - start + 1
            codon_start = start + (relative_pos - 1) // 3 * 3
            codon_end = codon_start + 3
            codon_pos = (relative_pos - 1) % 3

            ref_codon = str(ref_genome.seq[codon_start:codon_end]).upper()

            if len(ref_codon) == 3:  
                mut_codon = list(ref_codon)
                mut_codon[codon_pos] = major_base
                mut_codon = "".join(mut_codon)

                ref_aa = translate_codon(ref_codon)
                mut_aa = translate_codon(mut_codon)

                nonsynonymous = int(ref_aa != mut_aa)
                synonymous = int(ref_aa == mut_aa)
            else:
                # case where the codon is incomplete
                nonsynonymous = 0
                synonymous = 0
        else:
            nonsynonymous = 0
            synonymous = 0

        nonsynonymous_list.append(nonsynonymous)
        synonymous_list.append(synonymous)

    # Add the new columns to the filtered variants DataFrame
    filtered_variants['nonsynonymous'] = nonsynonymous_list
    filtered_variants['synonymous'] = synonymous_list
    print("Annotation complete")
    return filtered_variants

  

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Calculate number of nonsynonymous and synonymous sites in coding regions.")
    parser.add_argument("--output_file", type=str, required=True, help="Name of output CSV file")
    parser.add_argument("--ref_file", type=str, required=True, help="Name of reference FASTA file")
    parser.add_argument("--genbank_file", type=str, required=True, help="Path to the GenBank file")
    parser.add_argument("--timo_file", type=str, required=True, help="Path to output timo file")
    parser.add_argument("--maf", type=float, help="Minor allele frequency threshold", default=0.2)
    args = parser.parse_args()

    coding_regions = parse_genbank(args.genbank_file)
    filtered_data = parse_timo_output(args.timo_file, ref_fasta=args.ref_file, minor_threshold=args.maf)
    annotated_variants = annotate_timo_output(filtered_data, ref_genome=args.ref_file, coding_regions=coding_regions)
    annotated_variants.to_csv(args.output_file, index=False, header=True)
  
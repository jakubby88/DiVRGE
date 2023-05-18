## Written by: Kate Johnson
"""
INPUT: reference fasta used for alignment
OUTPUT: GTF file formatted to work with DiVRGE specific use

to run:
python3 make_gtf.py \
        -r ~/Lab/DiVRGE/testing/reference/H9N2_HK_1999_CDS.fasta \
        -s ./ \
        -p H9N2
"""
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--ref','-r', required=True,help='path & filename for ref')
parser.add_argument('--savedir','-s', required=True,help='directory path to save output gtf')
parser.add_argument('--prefix','-p', required=True,help='what you want to name gtf')
args = parser.parse_args()

def read_fasta(fp):
    """
    INPUT: Fasta file reference
    """
    name, seq = None, []
    for line in fp:
        line = line.rstrip()
        if line.startswith(">"):
            if name: yield (name, ''.join(seq))
            name, seq = line, []
        else:
            seq.append(line)
    if name: yield (name, ''.join(seq))

def open_fasta(filename):
    """
    INPUT: fasta file reference
    OUTPUT: dictionary with segment name as key and sequence of segment as value
    """
    segdict = {}
    with open(filename) as fp:
        for name, seq in read_fasta(fp):
            segdict[name[1:]] = seq
    return segdict

ref = open_fasta(args.ref)
filename = "{0}/{1}.gtf".format(args.savedir, args.prefix)

print("")
print("Input reference: {0}".format(args.ref))
print("")

# seqname = chrm name
# source: name of program
# feature: feature type (gene, exon)
# start: start of feature (1 = first nt of feature)
# end: end of feature
# score: floating pt value (?)
# straind + forward, - reverse
# frame: 0: first base of codon, 1=second base of codon, 3=third base
# https://useast.ensembl.org/info/website/upload/gff.html

with open(filename, 'w') as gtf:
    for SEGMENT in ref:
        print(SEGMENT, len(ref[SEGMENT]))
        long_string = 'gene_id "{0}"; gene_name "{0}"; gene_biotype "protein_coding"'.format(SEGMENT)
        gtf.write(SEGMENT + "\t" + 'divrge' + "\t" + "exon" + "\t" +  "1" + "\t" + "{0}".format(len(ref[SEGMENT])) + "\t" + "." + "\t" + "+" + "\t" + "." + "\t" +  long_string + "\n")
        #gtf.write(SEGMENT + "\t" + 'divrge' + "\t" + "exon" + "\t" + "1" + "\t" + "{0}".format(len(ref[SEGMENT])) + "\t" + "." + "\t" + "+" + "\t" + "." + "\t" +"\n")

gtf.close()
print("")
print("finished, check {0} for gtf".format(args.savedir))

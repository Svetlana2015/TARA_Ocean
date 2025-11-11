#Script by Duncan Sussfeld for Metagenomics analysis protocol, year 2022-2023
#Contact at duncan.sussfeld@gmail.com if needed

import os
import argparse
import pandas as pd
from collections import Counter
from ete3 import NCBITaxa
ncbi = NCBITaxa()


# Parameter setup
modes = ["group-subset", "otu-rep", "join-taxonomy"]
parser = argparse.ArgumentParser()
parser.add_argument("mode", type=str, choices=modes)
parser.add_argument("parameters", type=str, nargs="*")


def group_subset(seq_list_file, groups_file):
    """
    INPUT seq_list_file: text file listing sequence names for which to extract group information (1 per line)
    INPUT groups_file: text file (Mothur group file format) with group information for all sequences
    OUTPUT: writes filtered group file for provided sequence list
    """

    # all sequences in seq_list_file
    seqs = {l.strip() for l in open(seq_list_file).readlines()}

    # parse groups_file to extract info for all sequences in seqs
    basename, ext = os.path.splitext(groups_file)
    with open(basename + "_filtered" + ext, "w") as out_f:
        with open(groups_file) as in_f:
            for line in in_f:
                seq = line.split()[0]
                if seq in seqs:
                    out_f.write(line)
                    seqs.remove(seq)

    if len(seqs) != 0:
        print("WARNING: could not find group for the following {} sequences:".format(len(seqs)))
        print(" ".join(seqs))

    print("Filtered group file written at", basename + "_filtered" + ext)


def otu_rep(clusters_file, fasta_file):
    """
    INPUT clusters_file: Mothur-generated .list file (only one cutoff value)
    INPUT fasta_file: sequence file to extract representative sequences from
    OUTPUT: "oturep.fa" fasta file with one representative sequence per OTU, format >OtuXXXX|<seq_name>
    """
    rep_seqs = {}
    df = pd.read_table(clusters_file).T
    df.drop(index=["label", "numOtus"], inplace=True)
    for otu in df.index:
        seqs = df.at[otu, 0].split(",")
        rep_seq = seqs[0]
        rep_seqs[rep_seq] = otu

    with open(fasta_file) as in_f:
        out_path = os.path.join(os.path.dirname(fasta_file), "oturep.fa")
        with open(out_path, "w") as out_f:
            for line in in_f:

                if line.startswith(">"):
                    write_mode = False
                    seq_name = line.split()[0][1:]
                    if seq_name in rep_seqs:
                        write_mode = True
                        out_f.write(">" + rep_seqs[seq_name] + "|" + seq_name + "\n")
                else:
                    if write_mode:
                        out_f.write(line)

    print("Representative sequences written at oturep.fa")


def join_taxonomy(aln_file, ref_tax_file):
    """
    INPUT aln_file: BLAST alignment file where first two columns are qseqid and sseqid
    INPUT ref_tax_file: tab-delimited taxonomy annotation file for reference sequences
    OUTPUT: "oturep_tax_ids.txt" list of tax ids matched by OTUs
    OUTPUT: "oturep_tax_abundance.txt" abundance of taxa matched by OTUs
    """

    aln = pd.read_table(aln_file, names="qseqid sseqid evalue pident cov".split())
    tax = pd.read_table(ref_tax_file, names=["sseqid", "tax_id"], index_col=0)

    df = aln.join(tax, on="sseqid")

    abundance = dict(Counter(df["tax_id"]))

    with open("oturep_tax_ids.txt", "w") as out_f:
        out_f.write(" ".join([str(tax_id) for tax_id in abundance.keys()]))

    with open("oturep_tax_abundance.txt", "w") as out_f:
        out_f.write("DATASET_SIMPLEBAR\n")
        out_f.write("SEPARATOR TAB\n")
        out_f.write("DATASET_LABEL\tAbundance\n")
        out_f.write("COLOR\t#cc0099\n")
        out_f.write("DATA\n")
        taxid2name = ncbi.get_taxid_translator(list(abundance.keys()))
        for tax_id, abund in abundance.items():
            tax_name = taxid2name[tax_id].replace("=", "_")
            out_f.write(tax_name + " - " + str(tax_id) + "\t" + str(abund) + "\n")



if __name__ == "__main__":

    args = parser.parse_args()

    if args.mode == "group-subset":
        if len(args.parameters) != 2:
            raise IndexError("Expected exactly two arguments for group-subset mode!")
        else:
            group_subset(args.parameters[0], args.parameters[1])

    elif args.mode == "otu-rep":
        if len(args.parameters) != 2:
            raise IndexError("Expected exactly two arguments for otu-rep mode!")
        else:
            otu_rep(args.parameters[0], args.parameters[1])

    elif args.mode == "join-taxonomy":
        if len(args.parameters) != 2:
            raise IndexError("Expected exactly two arguments for join-taxonomy mode!")
        else:
            join_taxonomy(args.parameters[0], args.parameters[1])

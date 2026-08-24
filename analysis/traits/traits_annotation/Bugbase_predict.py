#!/bin/python
"""
Annotate full-length 16S rRNA sequences from the Greengenes2 database with
BugBase phenotype predictions.

Background
----------
BugBase provides precalculated phenotypic traits (e.g. Gram staining, oxygen
tolerance, biofilm formation, pathogenicity, mobile element content) for the
16S sequences contained in its own reference database. To transfer these
annotations onto Greengenes2 sequences, each query sequence is mapped to a
BugBase reference sequence by global alignment with vsearch. Only sequences
with an exact match to a BugBase reference inherit its phenotype; sequences
without such a match remain unannotated and are simply dropped.

Inputs
------
data/feces_seq_16S_new.fasta        Query 16S sequences (Greengenes2).
data/vsearch.out                    vsearch global alignment output in
                                    blast6out format (12 columns).
data/default_traits_precalculated.txt
                                    BugBase trait table, indexed by BugBase
                                    reference sequence ID.

Procedure
---------
1. Read the query FASTA and record the full length of every sequence.
2. Load the vsearch hits and compute, for each hit, the query coverage
   (alignment_length / query length).
3. Keep only hits that represent a perfect, full-length match:
       - alignment starts at the first base of the query (q.start == 1)
       - coverage == 1 (the entire query is aligned)
       - identity == 100%
   Queries failing these criteria are considered unannotatable and are
   excluded from the output.
4. Extract the BugBase traits of the retained reference sequences and relabel
   the rows with the corresponding query IDs.

Output
------
data/traits_precalculated.txt       Tab-separated BugBase trait table indexed
                                    by query (Greengenes2) sequence ID, ready
                                    to be used as the precalculated-trait file
                                    in downstream BugBase analysis.

Note
----
Because the filtering requires 100% identity over the full query length, a
single query may in principle match several BugBase references; likewise, only
the subset of Greengenes2 sequences present in BugBase will appear in the
output.
"""
import pandas as pd
import numpy as np
import re

from Bio import SeqIO

def get_sequence_lengths(fasta_file):
    seq_lengths = {}
    for record in SeqIO.parse(fasta_file, "fasta"):
        seq_lengths[record.id] = len(record.seq) 
    return seq_lengths

fasta_path = "../data/feces_seq_16S_silva.fasta"
length_dict = get_sequence_lengths(fasta_path)

colnames = ["query_id", "refer_id", "identity", "alignment_length", "mismatches", "gap_openings", "q.start",
            "q.end", "s.start", "s.end", "e-value", "bit_score"]
vsearch_out = pd.read_csv("../data/vsearch.out", sep="\t", header=None)
vsearch_out.columns = colnames
vsearch_out.loc[:, "length"] = [length_dict[i] for i in vsearch_out.query_id.values]
vsearch_out.loc[:, "coverage"] = vsearch_out.alignment_length.values / vsearch_out.length.values
vsearch_out = vsearch_out.loc[vsearch_out["q.start"] == 1]
vsearch_out = vsearch_out.loc[vsearch_out["coverage"] == 1]
vsearch_out = vsearch_out.loc[vsearch_out.identity == 100]

fid = vsearch_out.refer_id.values
file_path = "Bugbase_database.txt"
bugbase_data = pd.read_csv(file_path, sep = "\t", index_col = 0)
bugbase_data = bugbase_data.loc[fid]
bugbase_data.index = vsearch_out.query_id.values
bugbase_data.to_csv("traits_predict_Bugbase.txt", sep="\t")

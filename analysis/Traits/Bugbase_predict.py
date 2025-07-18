#!/bin/python

import pandas as pd
import numpy as np
import re

from Bio import SeqIO

def get_sequence_lengths(fasta_file):
    seq_lengths = {}
    for record in SeqIO.parse(fasta_file, "fasta"):
        seq_lengths[record.id] = len(record.seq)  # 存储序列ID和长度
    return seq_lengths

fasta_path = "Data/feces_seq_16S_new.fasta"
length_dict = get_sequence_lengths(fasta_path)

colnames = ["query_id", "refer_id", "identity", "alignment_length", "mismatches", "gap_openings", "q.start",
            "q.end", "s.start", "s.end", "e-value", "bit_score"]
vsearch_out = pd.read_csv("Data/vsearch.out", sep="\t", header=None)
vsearch_out.columns = colnames
vsearch_out.loc[:, "length"] = [length_dict[i] for i in vsearch_out.query_id.values]
vsearch_out.loc[:, "coverage"] = vsearch_out.alignment_length.values / vsearch_out.length.values
vsearch_out = vsearch_out.loc[vsearch_out["q.start"] == 1]
vsearch_out = vsearch_out.loc[vsearch_out["coverage"] == 1]
vsearch_out = vsearch_out.loc[vsearch_out.identity == 100]

fid = vsearch_out.refer_id.values
file_path = "Data/default_traits_precalculated.txt"
bugbase_data = pd.read_csv(file_path, sep = "\t", index_col = 0)
bugbase_data = bugbase_data.loc[fid]
bugbase_data.index = vsearch_out.query_id.values
bugbase_data.to_csv("Data/traits_precalculated.txt", sep="\t")

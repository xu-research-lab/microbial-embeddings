import sys
import biom
import subprocess
from Bio import SeqIO

# input
biom_file = sys.argv[1]
input_file = "/beegfs/db/greengenes/gg_13_8_otus/rep_set/99_otus.fasta"
fasta_file = sys.argv[2]
temp_path = sys.argv[3]
tree_path = sys.argv[4]

# 定义需要选取的序列名称
# table = biom.load_table(biom_file)
# selected_names = table.ids(axis='observation')
# selected_sequences = []
# for record in SeqIO.parse(input_file, "fasta"):
#     if record.id in selected_names:
#         selected_sequences.append(record)
# # 保存选取的序列为fasta文件
# SeqIO.write(selected_sequences, fasta_file, "fasta")
# 
# ## generate the tree by qiime2
# subprocess.call(f"qiime tools import \
#                  --type 'FeatureData[Sequence]' \
#                  --input-path {fasta_file} \
#                  --output-path {temp_path}/rep-seqs.qza", shell=True)

subprocess.call(f"qiime phylogeny align-to-tree-mafft-fasttree \
                 --i-sequences {temp_path}/rep-seqs.qza \
                 --o-alignment {temp_path}/aligned-rep-seqs.qza \
                 --o-masked-alignment {temp_path}/masked-aligned-rep-seqs.qza \
                 --o-tree {temp_path}/unrooted-tree.qza \
                 --o-rooted-tree {temp_path}/rooted-tree.qza", shell=True)

subprocess.call(f"qiime tools export \
                  --input-path {temp_path}/rooted-tree.qza \
                  --output-path {tree_path}", shell=True)

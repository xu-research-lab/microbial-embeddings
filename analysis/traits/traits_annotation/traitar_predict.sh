#!/bin/bash  

#SBATCH --nodes=1  
#SABTCH --ntasks=1  
#SBATCH --cpus-per-task=28  
#SBATCH -p cu  
#SBATCH --array=1%1 
#SBATCH --mem 250G  
#SBATCH --exclude=cu01,cu03  

module load miniconda/4.9.2  

source activate mpa

pfam_dir=/beegfs/db/pfam-33.1/
input_dir="genome_dir"
sample2file=samples.txt
output_dir="outpt_dir"
traitar phenotype ${pfam_dir} ${input_dir} ${sample2file} from_nucleotides ${output_dir} -o -c 28 -x 14

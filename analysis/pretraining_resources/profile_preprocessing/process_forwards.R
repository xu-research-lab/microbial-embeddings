# With many thanks to https://astrobiomike.github.io/amplicon/dada2_workflow_ex

library(dada2)
library(stringr)
packageVersion("dada2")

trimLeft = 0
trimRight = 0
args = commandArgs(trailingOnly=TRUE)
setwd(args[1])
if (length(args) > 1) {
    if (args[2] == 'trim') {
        trimLeft = 24
    }
}
log <- function(message) print(paste(date(), message))

# first we're setting a few variables we're going to use
# one with all sample names, by scanning our "samples" file we made earlier
samples_file <- scan("rum_sample.txt", what="character")

has_underscore_1 <- sapply(samples_file, function(filename) str_detect(filename, fixed("_1")))  
paried_samples <- samples_file[has_underscore_1]
single_samples <- samples_file[sapply(samples_file, function(filename) !grepl("_[12]", filename))]

# one holding the file names of all the forward reads
samples_file <- c(single_samples, paried_samples)
forward_reads <- paste0("rawdata/", samples_file)

sample_paired <- sub("_1.*", "", paried_samples)
sample_single <- sapply(strsplit(single_samples, "\\.fastq"), function(x) x[1])
samples <- c(sample_single, sample_paired)
# 检查哪些元素是重复的（TRUE = 重复，FALSE = 首次出现）
is_duplicated <- !duplicated(samples)
samples <- samples[is_duplicated]
forward_reads <- forward_reads[is_duplicated]
filtered_forward_reads <- paste0("intermediate/", samples, ".R1.filtered.fastq.gz")


#########################
# Quality filtering
#########################

# Filter for quality. Most important because it does a few
# DADA2-specific things (throwing out 'N' bases, for example)
# that become important later.
log('Filtering...')

filtered_out <- filterAndTrim(forward_reads, filtered_forward_reads,
                              truncQ=0, rm.phix=TRUE, multithread=TRUE,
                              verbose=TRUE, trimLeft=trimLeft, trimRight=trimRight)

# Then we limit the list of filtered fastq files to include
# only the ones that actually had reads pass the filter:
filtered_forward_reads <- filtered_forward_reads[file.exists(filtered_forward_reads)]
LOST_READS <- length(filtered_forward_reads) < length(forward_reads)

# revise the list of samples to only include those
# that actually have reads now:
samples <- gsub("^intermediate/", "", filtered_forward_reads)  # 去掉开头
samples <- gsub("\\.R1\\.filtered\\.fastq\\.gz$", "", samples)  # 去掉结尾

#########################
# Building error models
#########################
log('Building forward error model...')
err_forward_reads <- learnErrors(filtered_forward_reads, multithread=TRUE)

log('Plotting error models...')
pdf('temp/forward_error_model.pdf')
plotErrors(err_forward_reads, nominalQ=TRUE)
dev.off()

saveRDS(err_forward_reads, 'temp/err_forward_reads.rds')
# err_forward_reads <- readRDS('../temp/err_forward_reads.rds')

#########################
# Generate count table
#########################
ddF <- vector("list", length(samples))
names(ddF) <- samples

for(sam in samples) {
  cat("Processing:", sam, "\n")
  derepF <- derepFastq(paste("intermediate/", sam, ".R1.filtered.fastq.gz", sep=""))
  ddF[[sam]] <- dada(derepF, err=err_forward_reads, multithread=28)
}
rm(derepF)

seqtab <- makeSequenceTable(ddF)

# Get rid of really short sequences that can't practically be used
# to assign taxonomy:
seqtab.noshort <- seqtab[,nchar(colnames(seqtab)) > 49]
diff <- length(colnames(seqtab)) - length(colnames(seqtab.noshort))
log(paste('Removed',diff,'ASVs for being too short.'))

# check for chimeras
log('Removing bimeras...')
seqtab.nochim <- removeBimeraDenovo(seqtab.noshort, verbose=TRUE)

#########################
# Check reads dropped at each step
#########################
log('Writing ASV table...')
write.table(seqtab.nochim, "results/ASV.tsv",
            sep="\t", quote=F, col.names=NA)
saveRDS(seqtab.nochim, 'asv.rds')
log('ASVs recorded.')

# log('Assigning taxonomy...')
# taxa <- assignTaxonomy(seqtab.nochim, "/code/silva_nr99_v138_train_set.fa.gz", multithread=8, tryRC=T)

asv_seqs <- colnames(seqtab.nochim)
asv_headers <- vector(dim(seqtab.nochim)[2], mode="character")
for (i in 1:dim(seqtab.nochim)[2]) {
  asv_headers[i] <- paste(">ASV", i, sep="_")
}
# extract fasta:
asv_fasta <- c(rbind(asv_headers, asv_seqs))
# count table:
asv_tab <- t(seqtab.nochim)
row.names(asv_tab) <- sub(">", "", asv_headers)
# tax table:
# asv_tax <- taxa
# row.names(asv_tax) <- sub(">", "", asv_headers)


## output
log('Writing output files...')
write(asv_fasta, "results/ASVs.fa")
write.table(asv_tab, "results/ASVs_counts.tsv",
            sep="\t", quote=F, col.names=NA)
# write.table(asv_tax, "results/ASVs_taxonomy.tsv",
#             sep="\t", quote=F, col.names=NA)

getN <- function(x) sum(getUniques(x))

print('Calculating summary stats...')
# making a little table
forwd_val <- sapply(ddF, getN)
nochim_val <- rowSums(seqtab.nochim)
length_val <- rowSums(seqtab.noshort)
chim_removed_val <- round(((length_val-nochim_val)/forwd_val)*100, 1)

if(LOST_READS) {
  # if some samples had zero reads go through the filter,
  # we can't display the total input reads, because that
  # column has more entries than the rest of the columns,
  # which exclude the filtered-out samples
  summary_tab <- data.frame(forwd=forwd_val,
                          length=length_val,
                          nonchim=nochim_val,
                          chim_perc=chim_removed_val,
                          retained_perc=round((100*nochim_val)/forwd_val, 1))
} else {
  summary_tab <- data.frame(dinput=filtered_out[,1],
                          filter=filtered_out[,2], forwd=forwd_val,
                          length=length_val,
                          nonchim=nochim_val,
                          chim_perc=chim_removed_val,
                          retained_perc=round((100*nochim_val)/filtered_out[,1], 1))
}

# OUTPUT
log('Writing summary output...')
write.table(summary_tab, "results/summary.tsv",
            sep="\t", quote=F, col.names=NA)

log('DONE!!!')


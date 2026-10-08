# https://bioconductor.statistik.tu-dortmund.de/packages/3.7/bioc/vignettes/ropls/inst/doc/ropls-vignette.pdf
## Standard PLS-DA (ropls) with a label-permutation test.
## The phylogenetically corrected variant is kept in run_plsda_phylo.Rd.
library(ropls)
library(tidyverse)
library(doParallel)
library(doRNG)

N_CORES <- 10

emb <- read.csv("../../data/social_niche_embedding_100.txt",
                row.names = 1, sep = " ", header = FALSE)

## Fit one trait; returns NULL when there are too few labelled taxa.
## A slim list is stored so the notebook can plot without refitting.
plain_plsda <- function(x_all, trait_vec, n_perm = 999, min_n = 20) {
    ids <- intersect(names(trait_vec)[!is.na(trait_vec)], rownames(x_all))
    if (length(ids) < min_n || length(unique(trait_vec[ids])) < 2) return(NULL)

    y <- factor(trait_vec[ids])
    m <- opls(as.matrix(x_all[ids, ]), y, predI = 2, crossvalI = 5, permI = n_perm,
              fig.pdfC = "none", info.txtC = "none")

    ## permMN row 1 is the model with the observed labels. Standard
    ## permutation p value (1 + b) / (1 + m) (Phipson & Smyth 2010);
    ## ropls' own pQ2 divides by m instead.
    perm <- m@suppLs$permMN
    p_perm <- function(col) (1 + sum(perm[-1, col] >= perm[1, col])) / nrow(perm)

    list(q2_obs  = m@summaryDF[, "Q2(cum)"],
         r2y_obs = m@summaryDF[, "R2Y(cum)"],
         pQ2     = p_perm("Q2(cum)"),
         pR2Y    = p_perm("R2Y(cum)"),
         scores  = as.data.frame(m@scoreMN),            # columns p1, p2
         n_taxa  = length(ids),
         class_n = table(y))
}

## Fit every trait of one source in parallel; each worker writes its own
## data/<source>/plsda_<trait>.rds, errors are reported by the master.
run_source <- function(traits_df, traits_name, outdir) {
    outdir <- normalizePath(outdir)
    err <- foreach(i = traits_name, .packages = "ropls",
                   .export = c("plain_plsda", "emb")) %dopar% {
        tryCatch({
            res <- plain_plsda(emb, setNames(as.character(traits_df[, i]), rownames(traits_df)))
            if (!is.null(res)) saveRDS(res, file.path(outdir, paste0("plsda_", i, ".rds")))
            NULL
        }, error = function(e) conditionMessage(e))
    }
    for (k in which(!vapply(err, is.null, logical(1))))
        message("FAILED ", outdir, " / ", traits_name[k], ": ", err[[k]])
}

cl <- makeCluster(N_CORES)
registerDoParallel(cl)
registerDoRNG(20240809)                 # reproducible permutations across workers


## ================================================================
## bugbase
## ================================================================
traits <- read.csv("data/traits_bugbase.csv", row.names = 1)
traits[traits == ""] <- NA
traits_name <- c("Oxygen_Preference", "Gram_Status")

run_source(traits, traits_name, "data/bugbase")


## ================================================================
## Traitar
## ================================================================
traits <- read.csv("data/trait_predcit.csv", row.names = 1)
traits[traits == 3] <- 1

### Aerobe Facultative Anaerobe
traits$Oxygen_Preference <- NA
traits[traits$Aerobe == 1, "Oxygen_Preference"] <- 1
traits[traits$Facultative == 1, "Oxygen_Preference"] <- 2
traits[traits$Anaerobe == 1, "Oxygen_Preference"] <- 3
traits[rowSums(traits[, c("Aerobe", "Facultative", "Anaerobe")]) != 1, "Oxygen_Preference"] <- NA

### Gram negative Gram positive
traits$Gram_Status <- NA
traits[traits$Gram.negative == 1, "Gram_Status"] <- 1
traits[traits$Gram.positive == 1, "Gram_Status"] <- 2
traits[rowSums(traits[, c("Gram.negative", "Gram.positive")]) != 1, "Gram_Status"] <- NA

### Cell shape
traits$cell_shape <- NA
traits[traits$Coccus == 1, "cell_shape"] <- 1
traits[traits$Bacillus.or.coccobacillus == 1, "cell_shape"] <- 2
traits[rowSums(traits[, c("Coccus", "Bacillus.or.coccobacillus")]) != 1, "cell_shape"] <- NA

traits_name <- setdiff(colnames(traits),
                       c("Aerobe", "Facultative", "Anaerobe",
                         "Gram.negative", "Gram.positive",
                         "Coccus", "Bacillus.or.coccobacillus"))

run_source(traits, traits_name, "data/traitar")


## ================================================================
## bacDive
## ================================================================
traits <- read.csv("data/bacDive.csv")
traits <- traits[!duplicated(traits$X16s_ID), ]
rownames(traits) <- traits$X16s_ID

## bacDive is keyed by accession, the embedding by accession.start.end
df <- data.frame(accessions = str_split(rownames(emb), "\\.", simplify = TRUE)[, 1],
                 embed_id   = rownames(emb)) %>%
    filter(accessions %in% rownames(traits))
traits <- traits[df$accessions, ]
rownames(traits) <- df$embed_id

traits <- traits %>%
    mutate(Oxygen.Preference = case_when(anaerobe == 1               ~ "Anaerobic",
                                         facultative.anaerobe == 1   ~ "Facultatively",
                                         aerobe == 1                 ~ "Aerobic"))

## same missing values as load_bacdive() in traits_predict.ipynb
for (v in c("NA", "", "mixed", "variable", "+;NA", "no;yes", "negative;positive", "negative;variable"))
    traits[traits == v] <- NA
traits[traits == "-" | traits == 0] <- "no"
traits[traits == "+" | traits == 1] <- "yes"
traits[traits == "coccus-shaped"] <- "coccus"
traits[traits == "rod-shaped"]    <- "bacillus"
traits <- traits %>%
    mutate(cell_shape = if_else(cell_shape %in% c("coccus", "bacillus"), cell_shape, NA))

agg_bac <- read.csv("data/agg_bac.csv") %>%
    filter(level_2 %in% c("assimilation", "builds_acid_from"),
           terms %in% colnames(traits)) %>%
    distinct(terms)
traits <- traits[, c("gram_stain", "Oxygen.Preference", "cell_shape",
                     "spore_formation", "motility", agg_bac$terms)]
traits_name <- colnames(traits)[colSums(is.na(traits)) < nrow(traits)]

run_source(traits, traits_name, "data/bacdive")

stopCluster(cl)

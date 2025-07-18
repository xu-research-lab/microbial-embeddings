
# https://bioconductor.statistik.tu-dortmund.de/packages/3.7/bioc/vignettes/ropls/inst/doc/ropls-vignette.pdf
# ref: An approach for metabonomics data analysis applied on the plasma of Ginger water extract administered reserpine induced spleen deficiency rats
library(ropls)
library(tidyverse)
library(doParallel)
library(foreach)

### bugbase
traits <- read.csv("Data/traits_bugbase.csv", 
                   row.names = 1)
fid <- rownames(traits)
co_embedding <- read.csv("../../data/social_niche_embedding_100.txt",
                         row.names = 1, sep=" ", header = FALSE)

traits_name <- c("Oxygen_Preference", "Gram_Status")
for (i in traits_name){
    df =  co_embedding
    df$Trait = traits[, i]
    df <- df %>% filter(Trait != "")
    y <- as.factor(df$Trait)
    x <- as.matrix(co_embedding[rownames(df), ])
    plsda_model <- opls(x, y, predI=2, permI=999, crossvalI=5, fig.pdfC = 'none')
    saveRDS(plsda_model, file = paste0("Data/bugbase/model_", i,".rds"))
}

### Traitar
traits <- read.csv("Data/trait_predcit.csv", row.names = 1)
traits[traits == 3] <- 1

### Aerobe Facultative Anaerobe
traits$Oxygen_Preference <- rep(0, nrow(traits))
traits[traits$Aerobe == 1, "Oxygen_Preference"] <- 1
traits[traits$Facultative == 1, "Oxygen_Preference"] <- 2
traits[traits$Anaerobe == 1, "Oxygen_Preference"] <- 3
traits[rowSums(traits[, c("Aerobe", "Facultative", "Anaerobe")]) != 1, "Oxygen_Preference"] = NA

### Gram negative Gram positive
traits$`Gram_Status` <- rep(0, nrow(traits))
traits[traits$`Gram.negative` == 1, "Gram_Status"] <- 1
traits[traits$`Gram.positive` == 1, "Gram_Status"] <- 2
traits[rowSums(traits[, c("Gram.negative", "Gram.positive")]) != 1, "Gram_Status"] = NA

### Cell shape
traits$`cell_shape` <- rep(0, nrow(traits))
traits[traits$Coccus == 1, "cell_shape"] <- 1
traits[traits$`Bacillus.or.coccobacillus` == 1, "cell_shape"] <- 2
traits[rowSums(traits[, c("Coccus", "Bacillus.or.coccobacillus")]) != 1, "cell_shape"] = NA
fid <- rownames(traits)
co_embedding <- read.csv("../../data/social_niche_embedding_100.txt",
                         row.names = 1, sep=" ", header = FALSE)
co_embedding <- co_embedding[fid, ]

co_embedding_shuffled <- co_embedding %>% 
    mutate(across(everything(), ~ sample(.x)))

traits_name <- colnames(traits)
traits_name <- traits_name[!(traits_name %in% c("Aerobe", "Facultative", "Anaerobe",
                                                "Gram.negative", "Gram.positive",
                                                "Coccus", "Bacillus.or.coccobacillus"))]

n_cores <- detectCores() - 2  
cl <- makeCluster(n_cores)
registerDoParallel(cl)
output_dir <- "Data/traitar/"
foreach(i = traits_name, .packages = "ropls", .combine = 'c') %dopar% {
    df <- co_embedding
    df$Trait <- traits[, i]
    df <- df[!is.na(df$Trait), ]
    y <- as.factor(df$Trait)
    x <- as.matrix(co_embedding[rownames(df), ])
    plsda_model <- opls(x, y, predI = 2, crossvalI = 5, permI = 999, fig.pdfC = 'none')
    model_file <- paste0(output_dir, "model_", i, ".rds")
    saveRDS(plsda_model, file = model_file)
}

### bacDive
traits <- read.csv("Data/bacDive.csv")
traits <- traits[!duplicated(traits$`X16s_ID`), ]
rownames(traits) <- traits$`X16s_ID`

co_embedding <- read.csv("../../data/social_niche_embedding_100.txt", row.names = 1, sep=" ", header = FALSE)
accessions_num <- str_split(rownames(co_embedding), "\\.", simplify = TRUE)
accessions_num <- accessions_num[,1]
df <- data.frame(accessions=accessions_num, embed_id=rownames(co_embedding))
inter_id <- intersect(rownames(traits), df$accessions)
df <- df %>% filter(accessions %in% inter_id)
traits <- traits[df$accessions, ]
rownames(traits) <- df$embed_id
traits[traits == "NA"] <- NA
traits[traits == "-"] <- "no"
traits[traits == "+"] <- "yes"
traits[traits == ""] = NA
traits[traits == "+;NA"] = NA

traits$`Oxygen.Preference` <- rep(NA, nrow(traits))
traits <- traits %>% mutate(Oxygen.Preference = case_when(aerobe == 1 ~ "Aerobic", TRUE ~ Oxygen.Preference))
traits <- traits %>% mutate(Oxygen.Preference = case_when(`facultative.anaerobe` == 1 ~ "Facultatively", TRUE ~ Oxygen.Preference))
traits <- traits %>% mutate(Oxygen.Preference = case_when(anaerobe == 1 ~ "Anaerobic", TRUE ~ Oxygen.Preference))

traits[traits == "mixed"] = NA
traits[traits == ""] = NA
traits[traits == "filament-shaped"] = NA
traits[traits == "oval-shaped"] = NA
traits[traits == "ovoid-shaped"] = NA
traits[traits == "sphere-shaped"] = NA
traits[traits == "spiral-shaped"] = NA

remove_id <- c("X16s_ID", "aerobe", "facultative.anaerobe", "anaerobe")
traits <- traits[, !(colnames(traits) %in% remove_id)]
traits <- traits[, colSums(is.na(traits)) < nrow(traits)]

fid <- rownames(traits)
co_embedding <- co_embedding[fid, ]

co_embedding_shuffled <- co_embedding %>% 
    mutate(across(everything(), ~ sample(.x)))

traits_name <- colnames(traits)
traits_name <- traits_name[10:length(traits_name)]
output_dir <- "Data/bacdive/"

# 并行主循环
foreach(i = traits_name, .packages = "ropls", .combine = c, .errorhandling = "pass") %dopar% {
  result <- list(success = FALSE, trait = i)
  tryCatch({
    df <- co_embedding
    df$Trait <- traits[, i]
    df <- df[!is.na(df$Trait), ]
    
    if ((length(unique(df$Trait)) > 1) & (nrow(df) > 20)) {
      y <- as.factor(df$Trait)
      x <- as.matrix(co_embedding[rownames(df), ])
      plsda_model <- opls(x, y, predI = 2, crossvalI = 5, permI = 999, fig.pdfC = 'none')
      saveRDS(plsda_model, file = paste0(output_dir, "model_", i, ".rds"))
      result$success <- TRUE
    }
  }, error = function(e) {
    message(sprintf("Error processing %s: %s", i, e$message))
  })
  return(result)
}

stopCluster(cl)

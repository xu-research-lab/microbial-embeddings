library(tidyverse)
library(pROC)  
library(glmmTMB)
library(reshape2)
library(cowplot)
library(patchwork)

### caculate p-value
get_p <- function(sne_values, abund_values){
    data_long <- data.frame(
        values = c(sne_values, abund_values),
        method = factor(rep(c("SNE", "Abund"), each = length(sne_values))),
        study_id = factor(rep(1:length(sne_values), times = 2))
    )
    
    ### values must be 0 < y < 1 in glmmTMB model
    n <- nrow(data_long)
    data_long$values <- (data_long$values * (n - 1) + 0.5) / n
    
    model <- glmmTMB(values ~ method + (1 | study_id), 
                     data = data_long, 
                     family = beta_family(link = "logit"))
    return(summary(model)$coefficients)
}

prob_df <- read.csv("xxxx/all_disease_study_predictions.csv")

study_crc <- c("PRJNA318004", "PRJDB11845", "PRJNA290926", "PRJNA824020", 
               "PRJEB6070", "PRJNA430990", "PRJEB36789")
study_ibd <- c("qiita_2538", "PRJNA324147", "qiita_1629", "PRJNA450340", 
               "RISK_PRISM_f", "PRJNA368966", "PRJNA422193", "PRJNA431126")

### fix FPR = 0.2
disease_type <- c("CRC", "IBD")
plot_list <- list()
n <- 1
for (disease in disease_type){
    if (disease == "CRC"){
        study_id <- study_crc
    } else {
        study_id <- study_ibd
    }
    fpr = 0.2
    
    tpr_sne <- c()
    tpr_rf <- c()
    for (i in study_id){
        df_pick <- prob_df %>% filter(Disease == disease) %>% filter(Study == i)
        
        ## caculate TPR for sne
        roc_obj <- roc(df_pick$TrueLabel, df_pick$SNE_Probability)
        interp_function <- approxfun(x = 1 - roc_obj$specificities, y = roc_obj$sensitivities, 
                                     method = "linear",
                                     rule = 2) 
        tpr_sne <- c(tpr_sne, interp_function(fpr))
        
        ## caculate TPR for rf
        roc_obj <- roc(df_pick$TrueLabel, df_pick$Abund_Probability)
        interp_function <- approxfun(x = 1 - roc_obj$specificities, y = roc_obj$sensitivities, 
                                     method = "linear",  
                                     rule = 2) 
        tpr_rf <- c(tpr_rf, interp_function(fpr))
    }
    df <- data.frame(tpr_sne=tpr_sne, tpr_rf=tpr_rf, group=disease, study=study_id)
    
    p <- get_p(df$tpr_sne, df$tpr_rf)$cond["methodSNE", "Pr(>|z|)"]
    
    if (p < 0.001){
        p <- "p < 0.001"
    } else {
        p <- paste0("p = ", round(p, 3))
    }
    df_long <- df %>% melt(id.vars=c("group", "study"))
    mean_sne <- mean(df$tpr_sne)
    mean_rf <- mean(df$tpr_rf)
    plot_list[[n]] <- ggplot() +
        geom_point(data=df_long, aes(x=study, y=value, shape=variable), size=3)+
        geom_segment(data = df, aes(x = study, xend = study,
                                        y = tpr_sne, yend = tpr_rf),
                     color = "gray60", linewidth = 1.5, alpha = 0.7) +
        scale_shape_manual(values = c(16, 21),  
                           labels = c("SNE", "RF"))+
        geom_hline( yintercept = mean_sne, color = "black", linetype = "solid") + # dashed
        geom_hline( yintercept = mean_rf, color = "black", linetype = "dashed") +
        theme_bw(base_size = 10) +
        labs(y = "TPR", x = "", shape = "", title=paste0("FPR ", fpr, ": ", p)) +  
        theme(axis.text.x = element_text(angle = 15, hjust = 1),
              text = element_text(face = "bold"),
              plot.title = element_text(hjust = 0.5),
              axis.title = element_text(face = "bold"),
              axis.text = element_text(face = "bold"))
    n <- n + 1
}

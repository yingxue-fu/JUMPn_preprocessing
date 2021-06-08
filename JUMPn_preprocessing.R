# load required pacakges
library(tidyverse)
library(gtools)
library(svMisc)

############################## prepare the data ################################
# read example dataset
tcell_proteomics <- read.csv(file = "example_data/Example_dataset_t_cell_proteomics.csv")
# extract the protein expression matrix
tcell_proteomics_mx <- as.matrix(tcell_proteomics[,4:11])
row.names(tcell_proteomics_mx) <- tcell_proteomics$Gene.Names
# draw boxplot to make sure the data are normalized among samples
boxplot(log2(tcell_proteomics_mx))

################################### DE analysis ################################
# create a empty data.frame 
DE_results <- data.frame(matrix(NA, nrow = dim(tcell_proteomics_mx)[1], ncol = 4))
colnames(DE_results) = c("lgfc_0h_2h", "lgfc_0h_8h", "lgfc_0h_16h", "aov_P_value")
# compute fold changes
## mean values
mean0h <- apply(tcell_proteomics_mx[,1:2], 1, mean)
mean2h <- apply(tcell_proteomics_mx[,3:4], 1, mean)
mean8h <- apply(tcell_proteomics_mx[,5:6], 1, mean)
mean16h <- apply(tcell_proteomics_mx[,7:8], 1, mean)
## fold changes
DE_results$lgfc_0h_2h <- foldchange(mean2h, mean0h) %>% 
                          foldchange2logratio(base=2) %>% round(2)
DE_results$lgfc_0h_8h <- foldchange(mean8h, mean0h) %>% 
                          foldchange2logratio(base=2) %>% round(2)
DE_results$lgfc_0h_16h <- foldchange(mean16h, mean0h) %>% 
                          foldchange2logratio(base=2) %>% round(2)
# calculate aov and extract P-values
n_proteins <- dim(tcell_proteomics_mx)[1]
for (i in 1:n_proteins) {
  # data
  p <- data.frame(expression_value = tcell_proteomics_mx[i,],
                  group = as.factor(rep(c("0h", "2h", "8h", "16h"), each = 2)))
  # compute the analysis of variance
  p_aov <- aov(expression_value ~ group, data = p)
  # exract P value
  DE_results$aov_P_value[i] <- summary(p_aov)[[1]][["Pr(>F)"]][1] %>% 
                                formatC(digits = 1, format = "e")
  # show progress
  progress(i, max.value = n_proteins)
  Sys.sleep(0.01)
  if (i == n_proteins) message("Done!")
}
# adjust P-values for multiple comparisons
DE_results$BH_adjusted_P_value <- DE_results$aov_P_value %>% 
                                  p.adjust(method = "BH") %>% 
                                  formatC(digits = 1, format = "e")
# define DE proteins (max absolute log2FC > 0.5 and P-adj < 0.01)
DE_protein_index <- (apply(abs(DE_results[,1:3]), 1, max) > 0.5) & 
                    (as.numeric(DE_results$BH_adjusted_P_value) < 0.01)
DE_protein_expression_matrix <- tcell_proteomics[DE_protein_index,]
# output the data used for JUMPn
write.csv(DE_protein_expression_matrix, 
          file = "output/DE_protein_expression_matrix.csv", row.names = F)

#################### coefficient of variation (CV) calculation #################
# CV calculation
prot_exp_cv <- apply(tcell_proteomics_mx, 1, function(x) sd(x) / mean(x))
# choose top 20% proteins with largest CV
cv_top20_index <- prot_exp_cv > quantile(prot_exp_cv, 0.8)
CVtop20_protein_expression_matrix <- tcell_proteomics[cv_top20_index,]
# output the data used for JUMPn
write.csv(CVtop20_protein_expression_matrix, 
          file = "output/CVtop20_protein_expression_matrix.csv", row.names = F)

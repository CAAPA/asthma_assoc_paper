pca <- read.csv("./data/input/ADPC_GWAS PCA_results_2.csv", head=T, stringsAsFactors = F)
pca$PATIENT <- gsub("1:", "", pca$PATIENT)
pca.pcs <- pca[,c(1,4:13,2),]
pca.phenos <- pca[,c(1,3)]
pca.phenos$COHORT <- NA
pca.phenos$COHORT[pca.phenos$POPN == "Barbados"] <- "BAGS"
pca.phenos$COHORT[pca.phenos$POPN == "CEU"] <- "CEU"
pca.phenos$COHORT[pca.phenos$POPN == "Chicago"] <- "CAG"
pca.phenos$COHORT[pca.phenos$POPN == "Detroit"] <- "SAPPHIRE"
pca.phenos$COHORT[pca.phenos$POPN == "JHS_ARIC"] <- "ARIC"
pca.phenos$COHORT[pca.phenos$POPN == "JHS_JHS"] <- "JHS"
pca.phenos$COHORT[pca.phenos$POPN == "JHU_AA_650Y"] <- "GRAAD"
pca.phenos$COHORT[pca.phenos$POPN == "JHU_AA_ABR"] <- "BRIDGE"
pca.phenos$COHORT[pca.phenos$POPN == "NIH"] <- "NIH"
pca.phenos$COHORT[pca.phenos$POPN == "UCSF_AA"] <- "SAGE"
pca.phenos$COHORT[pca.phenos$POPN == "UCSF_PR"] <- "GALA"
pca.phenos$COHORT[pca.phenos$POPN == "WakeForest"] <- "SARP"
pca.phenos$COHORT[pca.phenos$POPN == "YRI"] <- "YRI"
pca.phenos$FID <- pca.phenos$PATIENT
pca.phenos$IID <- pca.phenos$PATIENT
pca.phenos$ORDER <- 0
pca.phenos$ORDER[pca.phenos$COHORT == "CEU"] <- 1
pca.phenos$ORDER[is.na(pca.phenos$COHORT)] <- 2
pca.phenos$ORDER[pca.phenos$COHORT == "GRAAD"] <- 3
pca.phenos$ORDER[pca.phenos$COHORT == "BRIDGE"] <- 4
pca.phenos$ORDER[pca.phenos$COHORT == "CAG"] <- 5
pca.phenos$ORDER[pca.phenos$COHORT == "SAPPHIRE"] <- 6
pca.phenos$ORDER[pca.phenos$COHORT == "ARIC"] <- 7
pca.phenos$ORDER[pca.phenos$COHORT == "JHS"] <- 8
pca.phenos$ORDER[pca.phenos$COHORT == "SAGE"] <- 9
pca.phenos$ORDER[pca.phenos$COHORT == "NIH"] <- 10
pca.phenos$ORDER[pca.phenos$COHORT == "SARP"] <- 11
pca.phenos$ORDER[pca.phenos$COHORT == "BAGS"] <- 12
pca.phenos$ORDER[pca.phenos$COHORT == "GALA"] <- 13

write.table(pca.pcs[order(pca.phenos$ORDER),], "./data/input/genesis_pca.evec",  sep="\t", quote=F, row.names=F, col.names=F)
pca.phenos <- pca.phenos[order(pca.phenos$ORDER),c("FID", "IID", "COHORT")]
write.table(pca.phenos, "./data/input/genesis_pca.phe",  sep=" ", quote=F, row.names=F, col.names=F)

#Steps to run PCA on Genesis:
#- Change x and y axis lables (PC not PCA)
#- CEU should be yellow (YRI red, NA blue)
#- Source populations should be square, the rest circles

write.table(pca.phenos, "./data/input/admixture_results/Admixture.phe",  sep=" ", quote=F, row.names=F, col.names=F)

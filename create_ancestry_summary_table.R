#Open the output txt file and then paste into supplementary table 4

sites <- c("jhu_650y",
           "jhu_abr",
           "chicago",
           "jackson_aric",
           "jackson_jhs",
           "ucsf_sf",
           "washington", 
           "winston_salem",
           "jhu_bdos",
           "ucsf_pr")

cohorts <- c("GRAAD", 
             "BRIDGE",
             "CAG",
             "ARIC",
             "JHS",
             "SAGE II",
             "NIH",
             "SARP",
             "BAGS",
             "GALA II"
)

#First determine which column holds African ancestry
admixture.results <- read.table("./data/input/admixture_results/Admixture.3.Q", stringsAsFactors = F)
admixture.fam <- read.table("./data/input/admixture_results/Admixture.fam", stringsAsFactors = F)
pca.ids <- read.csv("./data/input/ADPC_GWAS PCA_results_2.csv", head=T, stringsAsFactors = F)[,c(1,3)]
yri.ids <- sub("1:", "", pca.ids$PATIENT[which(pca.ids$POPN == "YRI")])
yri.rows <- which(admixture.fam$V1 %in% yri.ids)
yri.col <- as.numeric(which(colMeans(admixture.results[yri.rows,]) > 0.8))

out.frame <- data.frame()
i <- 0
for (site in sites) {
  #Get the cohort name
  i <- i + 1
  cohort <- cohorts[i]
  
  #Get African ancestry IQR
  site.ids <- read.csv(paste0("./data/input/pca_by_cohorts/pca_", site, ".csv"),stringsAsFactors = F)[1]
  site.ids <- sub("1:", "", site.ids$PATIENT)
  site.ids <- site.ids[grep("^WG", site.ids)]
  site.rows <- which(admixture.fam$V1 %in% site.ids)
  iqr <- round(quantile(admixture.results[site.rows,yri.col])[2:4],2)
  iqr.str <- paste0(iqr[2], " [", iqr[1], " - ", iqr[3], "]")
  
  #Get proporion variance explained
  file.name <- paste0("./data/input/pca_by_cohorts/eigen_", site, ".csv")
  if (file.exists(file.name)) {
    eigen.vec <- read.csv(file.name)
    total <- sum(eigen.vec[,2])
    prop.var.1 <- round(eigen.vec[1,2]/total,4)
    prop.var.2 <- round(eigen.vec[2,2]/total,4)
    prop.var.3 <- round(eigen.vec[3,2]/total,4)
  } else {
    prop.var.1 <- 0
    prop.var.2 <- 0
    prop.var.3 <- 0
  }

  #Get PC association values
  file.name <- paste0("./data/input/pca_by_cohorts/pca_", site, ".csv")
  if (file.exists(file.name)) {
    pca.frame <- read.csv(file.name, stringsAsFactors = F)[,c(2,4:6)]
    pca.frame <- pca.frame[pca.frame$GROUP %in% c(1,2),]
    pca.frame$GROUP <- pca.frame$GROUP-1
    p.1 <- round(summary(glm(GROUP ~ PC1, family = binomial(link = "logit"), data=pca.frame))$coefficients[2,4],3)
    p.2 <- round(summary(glm(GROUP ~ PC2, family = binomial(link = "logit"), data=pca.frame))$coefficients[2,4],3)
    p.3 <- round(summary(glm(GROUP ~ PC3, family = binomial(link = "logit"), data=pca.frame))$coefficients[2,4],3)
  } else {
    p.1 <- 0
    p.2 <- 0
    p.3 <- 0
  }
  
  #Add to out frame
  out.frame <- rbind(out.frame, 
                     data.frame(cohort, iqr.str, prop.var.1, prop.var.2, prop.var.3, p.1, p.2, p.3))
}

write.table(out.frame, "./data/output/ancestry_distro_summary.txt",  sep="\t", quote=F, row.names=F, col.names=F)

#Get median african ancestry by ethnicity
bdos.ids <- read.csv("./data/input/pca_by_cohorts/pca_jhu_bdos.csv",stringsAsFactors = F)[1]
bdos.ids <- sub("1:", "", bdos.ids$PATIENT)
bdos.ids <- bdos.ids[grep("^WG", bdos.ids)]
bdos.rows <- which(admixture.fam$V1 %in% bdos.ids)
pr.ids <- read.csv("./data/input/pca_by_cohorts/pca_ucsf_pr.csv",stringsAsFactors = F)[1]
pr.ids <- sub("1:", "", pr.ids$PATIENT)
pr.ids <- pr.ids[grep("^WG", pr.ids)]
pr.rows <- which(admixture.fam$V1 %in% pr.ids)
aa.ids <- read.csv("./data/input/ADPC_GWAS PCA_results_2.csv",stringsAsFactors = F)[1]
aa.ids <- sub("1:", "", aa.ids$PATIENT)
aa.ids <- aa.ids[grep("^WG", aa.ids)]
aa.ids <- aa.ids[which(!(aa.ids %in% c(bdos.ids, pr.ids)))]
aa.rows <- which(admixture.fam$V1 %in% aa.ids)


median(admixture.results[aa.rows,yri.col])
median(admixture.results[bdos.rows,yri.col])
median(admixture.results[pr.rows,yri.col])

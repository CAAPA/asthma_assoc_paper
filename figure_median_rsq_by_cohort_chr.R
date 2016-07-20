library(ggplot2)

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
             "GALA II")

cohorts.frame <- data.frame(site=sites, cohort=cohorts)

draw.frame <- data.frame()
for (site in sites){
  frame <- read.table(paste0("/Volumes/Promise Pegasus/caapa_imputation_pipeline/data/output/", 
                             site, "/imputed_qc/rsq_summary.txt"), head=T, stringsAsFactors = F)
  frame$site <- site
  draw.frame <- rbind(draw.frame, frame)
}
draw.frame <- draw.frame[draw.frame$category %in% c("all", "all_lt", "all_gt"),]
draw.frame$category[draw.frame$category == "all"] <- "All SNPs"
draw.frame$category[draw.frame$category == "all_lt"] <- "SNPs MAF <= 0.005"
draw.frame$category[draw.frame$category == "all_gt"] <- "SNPs MAF > 0.005"
draw.frame <- draw.frame[draw.frame$chr != "1-22",]
draw.frame$chr <- factor(draw.frame$chr, levels=as.character(1:22))
draw.frame <- merge(draw.frame, cohorts.frame)
draw.frame$cohort <- factor(draw.frame$cohort, levels=cohorts)

pdf("./data/output/median_rsq/supp_fig_median_rsq.pdf", width = 14, height=7)
ggplot(data=draw.frame, aes(x=chr, y=rsq, group=cohort, colour=cohort)) +
  facet_wrap(~category) +
  geom_line() +
  geom_point() +
  labs(x="Chromosome") +
  labs(y="Rsq") +
  theme_bw() +
  theme(legend.title = element_blank()) 
dev.off()
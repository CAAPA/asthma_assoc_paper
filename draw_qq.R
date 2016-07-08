args <- commandArgs(trailingOnly = TRUE)

study.name <- args[1] #"BRIDGE"
in.file.name <- args[2] #"tmp_jhu_abr.txt"
out.file.name <- args[3] #"jhu_abr.png"
err.p.file.name <- args[4] #"jhu_small_p.txt"

plotQQ <- function(p.vals, study.name, plot.type) {
  observed <- sort(p.vals)
  lobs <- -(log10(observed))
  expected <- c(1:length(observed))
  lexp <- -(log10(expected / (length(expected)+1)))
  plot(c(0,10), c(0,10), col="red", lwd=1, type="l", 
       xlab="Expected (-logP)", ylab="Observed (-logP)", 
       xlim=c(0,10), ylim=c(0,10), 
       las=1, xaxs="i", yaxs="i", bty="l", main=paste0(study.name, " (" , plot.type, ")"))
  points(lexp, lobs, pch=23, cex=.1, bg="black")
  chisq1 <- qchisq(p.vals,1,lower.tail = F)
  inf.fact <- round(median(chisq1)/qchisq(0.5,1),4)
  text(7, 2, bquote(lambda == .(inf.fact)))
}

png(out.file.name, width=9.5, height=3.8, units = 'in', res = 200)
par(mfrow=c(1,3), cex=1.2)
results <- read.table(in.file.name, head=T, stringsAsFactors = F)
write.table(results[results$PVALUE < 1e-9,], err.p.file.name,
            sep="\t", quote=F, row.names=F, col.names=T)
p.vals <- results$PVALUE
freq <- results$MAF
plotQQ(p.vals, study.name, "All")
plotQQ(p.vals[freq >= 0.05], study.name, "Common")
plotQQ(p.vals[freq < 0.05], study.name, "Rare")
dev.off()

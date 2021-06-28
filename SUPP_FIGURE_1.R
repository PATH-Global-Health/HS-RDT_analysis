
rm(list = ls())
library(dplyr)
library(binom)

cols = adjustcolor(rev(c("deepskyblue2", "firebrick1")), 0.7)


# UGANDA 1
dd1 = read.csv("data/BIOME1_uganda.csv", stringsAsFactors = FALSE)

# MYANMAR 2
dd2 = read.csv("data/SMRU_Myanmar Compiled_csv.csv", stringsAsFactors = FALSE)


dd1$HRP2 = dd1$Quansys.HRP2.avg.conc..pg.mL.
dd1$HRP2[dd1$HRP2 == "<0.1"] = 0
dd1$HRP2[dd1$HRP2 == ">14600"] = 14600
dd1$HRP2 = as.numeric(dd1$HRP2)

dd1$pfPCR2 = rep(0, nrow(dd1))
dd1$pfPCR2[dd1$Pf.qPCR.designation %in% c("Pf")] = 1

dd1 = dd1[dd1$pfPCR2 == 1,]


dd2$HRP2 = dd2$HRP2_pg_ml
dd2$HRP2[which(dd2$HRP2 == "< 1.07")] = 0
dd2$HRP2[which(dd2$HRP2 == "> 16500.00")] = 16500
dd2$HRP2 = as.numeric(dd2$HRP2)
dd2 = dd2[!is.na(dd2$HRP2),]

dd2$qPCR= rep(0, nrow(dd2))
dd2$qPCR[which(dd2$qPCR.result %in% c("PF", "PF+PV"))] = 1

dd2 = dd2[dd2$qPCR == 1,]

dd1$rdt_ind = rep(NA, nrow(dd1))
dd1$rdt_ind[which(dd1$HS.RDT.pos.neg == "Pos" & dd1$RDT.pos.neg == "Pos")] = 1
dd1$rdt_ind[which(dd1$HS.RDT.pos.neg == "Pos" & dd1$RDT.pos.neg %in% c("Neg", "discordant"))] = 2
dd1$rdt_ind[which(dd1$HS.RDT.pos.neg %in% c("Neg", "") & dd1$RDT.pos.neg == "Pos")] = 3
dd1$rdt_ind[which(dd1$HS.RDT.pos.neg %in% c("Neg", "") & dd1$RDT.pos.neg %in% c("Neg", "discordant"))] = 4

dd1_a = dd1[dd1$Pf.qPCR.designation == "Pf",]
# dd1_a = dd1_a[-which(dd1_a$HS.RDT.pos.neg == ""),]
# dd1_a = dd1_a[-which(dd1_a$RDT.pos.neg == "discordant"),]
dd1_a = dd1_a[which(dd1_a$Quansys.HRP2.designation == "Pos"),]



dd2$rdt_ind = rep(NA, nrow(dd2))
# dd2$rdt_ind[which(dd2$Lab.HSRDT == "PF" & dd2$Lab.RDT %in% c("PF", "PF+PV"))] = 1
# dd2$rdt_ind[which(dd2$Lab.HSRDT == "PF" & dd2$Lab.RDT  %in% c("Neg", "PV"))] = 2
# dd2$rdt_ind[which(dd2$Lab.HSRDT == "Neg" & dd2$Lab.RDT %in% c("PF", "PF+PV"))] = 3
# dd2$rdt_ind[which(dd2$Lab.HSRDT == "Neg" & dd2$Lab.RDT %in% c("Neg", "PV"))] = 4

dd2$rdt_ind[which(dd2$Lab.HSRDT == "PF" & dd2$Lab.RDT %in% c("PF", "PF+PV"))] = 1
dd2$rdt_ind[which(dd2$Lab.HSRDT == "PF" & dd2$Lab.RDT  %in% c("Neg"))] = 2
dd2$rdt_ind[which(dd2$Lab.HSRDT == "Neg" & dd2$Lab.RDT %in% c("PF", "PF+PV"))] = 3
dd2$rdt_ind[which(dd2$Lab.HSRDT == "Neg" & dd2$Lab.RDT %in% c("Neg"))] = 4


dd2_a = dd2[!is.na(dd2$rdt_ind),]
dd2_a = dd2_a[dd2_a$qPCR == 1, ]
dd2_a = dd2_a[which(dd2_a$HRP2_pg_ml_pos == "positive"), ]

cols2 = c("darkorchid4", "firebrick1", "deepskyblue2", "grey80")


cx.axis = 0.7
cx.names = 0.8
cx.leg = 0.5
cx.main = 0.7
cx = 0.8
cx.lab = 0.7


tiff("plots/FINAL_submitted/SUPP_FIGURE_S1.tiff",width=170,height=82,units="mm",res=300, compression="lzw")

par(mfrow = c(1,2), xpd=FALSE, mar = c(4,4,2.3,1.5))

plot(log10(dd1_a$Pf.qRT.PCR..p.uL.), log10(dd1_a$HRP2), xaxt="n", yaxt="n", xlab=expression(paste("qRT-PCR parasite density (parasites/", mu, "L)")),
     ylab= "HRP2 concentration (pg/ml)", bg=adjustcolor(cols2[dd1_a$rdt_ind], 0.7), pch=21, ylim = c(0, log10(16500)), xlim=c(-2, 7),cex=cx, 
     main = "High transmission (Uganda)",
     cex.main = cx.main, cex.axis = cx.axis, cex.lab = cx.lab)
abline(v = seq(-2, 7), col="grey", lty=3)
axis(1, -2:7, 10^(-2:7), cex.axis = cx.axis)
axis(2, -2:6, 10^(-2:6), las=1, cex.axis = cx.axis)
legend("bottomright", c("U-RDT only", "co-RDT only", "U-RDT and co-RDT", "none"), pt.bg=adjustcolor(cols2[c(2,3,1,4)], 0.7), 
       pt.cex = 1, pch=21, bty = "n", cex = cx.leg)


plot(log10(dd2_a$Parasitemia.qPCR), log10(dd2_a$HRP2), xaxt="n", yaxt="n", xlab=expression(paste("qPCR parasite density (parasites/", mu, "L)")),
     ylab= "HRP2 concentration (pg/ml)", bg=adjustcolor(cols2[dd2_a$rdt_ind], 0.7), pch=21, ylim = c(0, log10(16500)),cex=cx, xlim=c(-2,7), 
     main = "Low transmission (Myanmar)",
     cex.main = cx.main, cex.axis = cx.axis, cex.lab = cx.lab)
abline(v = seq(-2, 7), col="grey", lty=3)
axis(1, -2:7, 10^(-2:7), cex.axis = cx.axis)
axis(2, -2:6, 10^(-2:6), las=1, cex.axis = cx.axis)

dev.off()

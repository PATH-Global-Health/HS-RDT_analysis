# analytical plots
.libPaths("C:/Users/hslater/Documents/R/win-library/4.0")
# read in myanmar and uganda data

rm(list = ls())
library(dplyr)
library(binom)
library(brms)

cols = adjustcolor(rev(c("deepskyblue2", "firebrick1")), 0.7)


# UGANDA 1
dd1 = read.csv("data/BIOME1_uganda.csv", stringsAsFactors = FALSE)

# MYANMAR 2
dd2 = read.csv("data/SMRU_Myanmar Compiled_csv.csv", stringsAsFactors = FALSE)


# clean data names and data types 

dd1$HRP2 = dd1$Quansys.HRP2.avg.conc..pg.mL.
dd1$HRP2[dd1$HRP2 == "<0.1"] = 0.05
dd1$HRP2[dd1$HRP2 == ">14600"] = 14600
dd1$HRP2 = as.numeric(dd1$HRP2)

dd1$pfPCR2 = rep(0, nrow(dd1))
dd1$pfPCR2[dd1$Pf.qPCR.designation %in% c("Pf")] = 1

dd1 = dd1[dd1$pfPCR2 == 1,]


dd2$HRP2 = dd2$HRP2_pg_ml
dd2$HRP2[which(dd2$HRP2 == "< 1.07")] = 0.05
dd2$HRP2[which(dd2$HRP2 == "> 16500.00")] = 16500
dd2$HRP2 = as.numeric(dd2$HRP2)
dd2 = dd2[!is.na(dd2$HRP2),]

dd2$qPCR= rep(0, nrow(dd2))
dd2$qPCR[which(dd2$qPCR.result %in% c("PF", "PF+PV"))] = 1

dd2 = dd2[dd2$qPCR == 1,]


# define thresholds

cut_seq = c(0, 1, 10, 100, 1000, 10000, 100000)
cut_seq2 = c(0, 1, 10, 100, 1000, 10000, 1000000)

dd1$hrp2_cut = cut(dd1$HRP2, cut_seq, include.lowest = TRUE, labels = 1:(length(cut_seq)-1))
dd2$hrp2_cut = cut(dd2$HRP2, cut_seq, include.lowest = TRUE, labels = 1:(length(cut_seq)-1))

# calculate number of positives and number of tests for each dataset
# for samples within each group

zz1 = dplyr::summarise(group_by(dd1, hrp2_cut),
                       hsrdt_pos = length(which(HS.RDT.pos.neg == "Pos")),
                       hsrdt_n = length(which(HS.RDT.pos.neg %in% c("Pos", "Neg", ""))),
                       rdt_pos = length(which(RDT.pos.neg == "Pos")),
                       rdt_n = length(which(RDT.pos.neg %in% c("Pos", "Neg", "discordant"))),
                       geom_mean_pcr = 10^(mean(log10(Pan.qRT.PCR..p.uL.+0.01))))

zz2 = dplyr::summarise(group_by(dd2, hrp2_cut),
                       hsrdt_pos = length(which(Lab.HSRDT == "PF")),
                       hsrdt_n = length(which(Lab.HSRDT %in% c("PF", "Neg"))),
                       rdt_pos = length(which(Lab.RDT == "PF")),
                       rdt_n = length(which(Lab.RDT %in% c("PF", "Neg"))),
                       geom_mean_pcr = 10^(mean(log10(Parasitemia.qPCR+0.01))))

# calculate binomial confidence intervals around all proportions

#uganda
bc_hsrdt = binom.confint(zz1$hsrdt_pos, zz1$hsrdt_n, method="exact")
bc_rdt = binom.confint(zz1$rdt_pos, zz1$rdt_n, method="exact")

#myanmar
bc_hsrdt2 = binom.confint(zz2$hsrdt_pos, zz2$hsrdt_n, method="exact")
bc_rdt2 = binom.confint(zz2$rdt_pos, zz2$rdt_n, method="exact")

# make nice names for plot legend
options(scipen=999)
nms = paste0(cut_seq[1:(length(cut_seq)-1)], "-",cut_seq[2:length(cut_seq)])


###################################################################################################################
###################################################################################################################

# set up binary variables for logistic regression

dd1$HS_RDT = if_else(dd1$HS.RDT.pos.neg == "Neg", 0,1)
dd1$co_RDT = if_else(dd1$RDT.pos.neg == "Pos", 1,0)

dd2 = dd2[dd2$Lab.HSRDT != "ND",]
dd2$HS_RDT = if_else(dd2$Lab.HSRDT == "PF", 1, 0)
dd2$co_RDT = if_else(dd2$Lab.RDT %in% c("PF","PF+PV"), 1, 0)

# run logisitic regression models

mod1 = brms::brm(
  HS_RDT ~ log10(HRP2),
  data = dd1,
  family = "bernoulli",
  prior = prior('normal(0,3)'),
  iter = 5000,
  chains = 1,
  cores = 4
)

mod2 = brms::brm(
  co_RDT ~ log10(HRP2),
  data = dd1,
  family = "bernoulli",
  prior = prior('normal(0,3)'),
  iter = 5000,
  chains = 1,
  cores = 4
)


mod3 = brms::brm(
  HS_RDT ~ log10(HRP2),
  data = dd2,
  family = "bernoulli",
  prior = prior('normal(0,3)'),
  iter = 5000,
  chains = 1,
  cores = 4
)

mod4 = brms::brm(
  co_RDT ~ log10(HRP2),
  data = dd2,
  family = "bernoulli",
  prior = prior('normal(0,3)'),
  iter = 5000,
  chains = 1,
  cores = 4
)
saveRDS(mod1, "fitted models/HS_RDT_uganda.RDS")
saveRDS(mod2, "fitted models/co_RDT_uganda.RDS")
saveRDS(mod3, "fitted models/HS_RDT_myanmar.RDS")
saveRDS(mod4, "fitted models/co_RDT_myanmar.RDS")

# read data in again and clean PCR data

# UGANDA 1
bb1 = read.csv("data/BIOME1_uganda.csv", stringsAsFactors = FALSE)

# MYANMAR 2
bb2 = read.csv("data/SMRU_Myanmar Compiled_csv.csv", stringsAsFactors = FALSE)


bb1$HRP2 = bb1$Quansys.HRP2.avg.conc..pg.mL.
bb1$HRP2[bb1$HRP2 == "<0.1"] = 0.05
bb1$HRP2[bb1$HRP2 == ">14600"] = 14600
bb1$HRP2 = as.numeric(bb1$HRP2)



bb2$HRP2 = bb2$HRP2_pg_ml
bb2$HRP2[which(bb2$HRP2 == "< 1.07")] = 0.05
bb2$HRP2[which(bb2$HRP2 == "> 16500.00")] = 16500
bb2$HRP2 = as.numeric(bb2$HRP2)
bb2 = bb2[!is.na(bb2$HRP2),]

# leave out P.spp and equivocal samples
bb2 = bb2[bb2$qPCR.result %in% c("Negative", "PF", "PF+PV"),]

bb1$HS_RDT = if_else(bb1$HS.RDT.pos.neg == "Neg", 0,1)
bb1$co_RDT = if_else(bb1$RDT.pos.neg == "Pos", 1,0)

bb2 = bb2[bb2$Lab.HSRDT != "ND",]
bb2$HS_RDT = if_else(bb2$Lab.HSRDT == "PF", 1, 0)
bb2$co_RDT = if_else(bb2$Lab.RDT %in% c("PF","PF+PV"), 1, 0)

bb2$Parasitemia.qPCR[is.na(bb2$Parasitemia.qPCR)] = 0

bb2$Parasitemia.qPCR = bb2$Parasitemia.qPCR/1e3

# run PCR logisitic regresison models

mod_pcr_1 = brms::brm(
  HS_RDT ~ log10(Pf.qRT.PCR..p.uL.+0.01),
  data = bb1,
  family = "bernoulli",
  prior = prior('normal(0,3)'),
  iter = 5000,
  chains = 1,
  cores = 4
)


mod_pcr_2 = brms::brm(
  co_RDT ~ log10(Pf.qRT.PCR..p.uL.+0.01),
  data = bb1,
  family = "bernoulli",
  prior = prior('normal(0,3)'),
  iter = 5000,
  chains = 1,
  cores = 4
)


mod_pcr_3 = brms::brm(
  HS_RDT ~ log10(Parasitemia.qPCR+0.0000002),
  data = bb2,
  family = "bernoulli",
  prior = prior('normal(0,3)'),
  iter = 5000,
  chains = 1,
  cores = 4
)


mod_pcr_4 = brms::brm(
  co_RDT ~ log10(Parasitemia.qPCR+0.0000002),
  data = bb2,
  family = "bernoulli",
  prior = prior('normal(0,3)'),
  iter = 5000,
  chains = 1,
  cores = 4
)

saveRDS(mod_pcr_1, "fitted models/HS_RDT_pcr_uganda.RDS")
saveRDS(mod_pcr_2, "fitted models/co_RDT_pcr_uganda.RDS")
saveRDS(mod_pcr_3, "fitted models/HS_RDT_pcr_myanmar.RDS")
saveRDS(mod_pcr_4, "fitted models/co_RDT_pcr_myanmar.RDS")


mod1 <- readRDS("fitted models/HS_RDT_uganda.RDS")
mod2 <- readRDS("fitted models/co_RDT_uganda.RDS")
mod3 <- readRDS("fitted models/HS_RDT_myanmar.RDS")
mod4 <- readRDS("fitted models/co_RDT_myanmar.RDS")



mod_pcr_1 <- readRDS("fitted models/HS_RDT_pcr_uganda.RDS")
mod_pcr_2 <- readRDS("fitted models/co_RDT_pcr_uganda.RDS")
mod_pcr_3 <- readRDS("fitted models/HS_RDT_pcr_myanmar.RDS")
mod_pcr_4 <- readRDS("fitted models/co_RDT_pcr_myanmar.RDS")


xseq2 = 10^seq(-1, 5, by=0.01)
xseq3 = 10^seq(-1, 6, by=0.01)



# make predictions

pred_ci_1 = fitted(mod1, newdata = data.frame(HRP2  = xseq2), summary = TRUE)
pred_ci_2 = fitted(mod2, newdata = data.frame(HRP2  = xseq2), summary = TRUE)

pred_ci_3 = fitted(mod3, newdata = data.frame(HRP2  = xseq2), summary = TRUE)
pred_ci_4 = fitted(mod4, newdata = data.frame(HRP2  = xseq2), summary = TRUE)


pred_ci_pcr_1 = fitted(mod_pcr_1, newdata = data.frame(Pf.qRT.PCR..p.uL.  = xseq3), summary = TRUE)
pred_ci_pcr_2 = fitted(mod_pcr_2, newdata = data.frame(Pf.qRT.PCR..p.uL.  = xseq3), summary = TRUE)

pred_ci_pcr_3 = fitted(mod_pcr_3, newdata = data.frame(Parasitemia.qPCR  = xseq3), summary = TRUE)
pred_ci_pcr_4 = fitted(mod_pcr_4, newdata = data.frame(Parasitemia.qPCR  = xseq3), summary = TRUE)




cx.axis = 0.8
cx.names = 0.8
cx.leg = 0.8
cx.main = 0.8
cx = 0.7

#tiff("plots/analytical_fig1_barplots3.tiff",width=270,height=240,units="mm",res=300, compression="lzw")
tiff("plots/FINAL_submitted/FIGURE_2.tiff",width=170,height=150,units="mm",res=300, compression="lzw")

par(mar=c(5.5,4,2.9,1.5), mfrow=c(2,2))
bp = barplot(rbind(bc_hsrdt$mean*100, bc_rdt$mean*100), beside=TRUE,col=cols, border=NA, names.arg = nms, 
             las=2, ylab="", 
             main = "High transmission (Uganda)", cex.names = cx.names, cex.axis = cx.axis,
             cex.main =cx.main)
segments(bp[1,], bc_hsrdt$lower*100, bp[1,], bc_hsrdt$upper*100, lwd=2, col="firebrick3")
segments(bp[2,], bc_rdt$lower*100, bp[2,], bc_rdt$upper*100, lwd=2, col="grey40")
mtext("HRP2 concentration (pg/ml)", 1, line=5, cex = cx)
mtext("Sensitivity (%)", 2, line=2.5, cex = cx)

legend("topleft", c("HS-RDT","co-RDT"), col = cols, pch=15, pt.cex = 1.6, bty="n", cex=cx.leg)
mtext("A", side = 3, adj = -0.18, line=1, cex=1, font=2)

bp = barplot(rbind(bc_hsrdt2$mean*100, bc_rdt2$mean*100), beside=TRUE,col=cols, border=NA, names.arg = nms,
             las=2, ylab="", 
             main = "Low transmission (Myanmar)", cex.names = cx.names, cex.axis = cx.axis,
             cex.main =cx.main)
segments(bp[1,], bc_hsrdt2$lower*100, bp[1,], bc_hsrdt2$upper*100, lwd=2, col="firebrick3")
segments(bp[2,], bc_rdt2$lower*100, bp[2,], bc_rdt2$upper*100, lwd=2, col="grey40")
mtext("HRP2 concentration (pg/ml)", 1, line=5, cex = cx)
legend("topleft", c("HS-RDT","co-RDT"), col = cols, pch=15, pt.cex = 1.6, bty="n", cex=cx.leg)
mtext("Sensitivity (%)", 2, line=2.5, cex = cx)
mtext("B", side = 3, adj = -0.18, line=1, cex=1, font=2)


cols2 = adjustcolor(rev(c("deepskyblue2", "firebrick1")), 0.7)

par(mar=c(4,4,3.3,1.5))

plot(log10(xseq2), pred_ci_1[,1], type="l", xaxt="n", las=1, xlab="",
     main = "Probability of positivity by HRP2 concentration",
     ylab = "", lwd=3, col=cols2[1],  cex.axis = cx.axis,
     cex.main =cx.main)
abline(h = seq(0,1, by=0.1), col="grey80", lty=3)
abline(v= seq(-1, 5, by=1), col="grey80", lty=3)
mtext("HRP2 concentration (pg/ml)", 1, line=3, cex = cx)
mtext("Probabilty of HS-RDT or co-RDT positivity", 2, line=3, cex = cx)

lines(log10(xseq2), pred_ci_1[,1], lwd=3, col=cols2[1])
polygon(c(log10(xseq2), rev(log10(xseq2))), c(pred_ci_1[,3], rev(pred_ci_1[,4])), 
        col=adjustcolor(cols2[1], 0.3), border=NA)

lines(log10(xseq2), pred_ci_2[,1], lwd=3, col=cols2[2])
polygon(c(log10(xseq2), rev(log10(xseq2))), c(pred_ci_2[,3], rev(pred_ci_2[,4])), 
        col=adjustcolor(cols2[2], 0.3), border=NA)

lines(log10(xseq2), pred_ci_3[,1], lwd=3, col=cols2[1], lty=5)
polygon(c(log10(xseq2), rev(log10(xseq2))), c(pred_ci_3[,3], rev(pred_ci_3[,4])), 
        col=adjustcolor(cols2[1], 0.3), border=NA)

lines(log10(xseq2), pred_ci_4[,1], lwd=3, col=cols2[2], lty=5)
polygon(c(log10(xseq2), rev(log10(xseq2))), c(pred_ci_4[,3], rev(pred_ci_4[,4])), 
        col=adjustcolor(cols2[2], 0.3), border=NA)

axis(1, -1:5, 10^(-1:5), cex.axis = cx.axis)

# legend("topleft", c("HS-RDT, Uganda", "HS-RDT, Myanmar", "co-RDT, Uganda", "co-RDT, Myanmar"),
#        col=rep(cols2, each = 2), lwd=3, lty=c(1,2,1,2), bty="n", cex = cx)

mtext("C", side = 3, adj = -0.18, line=1, cex=1, font=2)

##########################################################################################################################

plot(log10(xseq3), pred_ci_pcr_1[,1], type="l", xaxt="n", las=1, 
     xlab="", main = "Probability of positivity by parasite density",
     ylab = "", lwd=3, col=cols2[1], ylim=c(0,1),  cex.axis = cx.axis,
     cex.main =cx.main)
abline(h = seq(0,1, by=0.1), col="grey80", lty=3)
abline(v= seq(-1, 8, by=1), col="grey80", lty=3)

mtext(expression(paste("Parasite density (parasites per ",mu,"l)")), 1, line=3, cex = cx)
mtext("Probabilty of HS-RDT or co-RDT positivity", 2, line=3, cex = cx)

lines(log10(xseq3), pred_ci_pcr_1[,1], lwd=3, col=cols2[1])
polygon(c(log10(xseq3), rev(log10(xseq3))), c(pred_ci_pcr_1[,3], rev(pred_ci_pcr_1[,4])), 
        col=adjustcolor(cols2[1], 0.3), border=NA)

lines(log10(xseq3), pred_ci_pcr_2[,1], lwd=3, col=cols2[2])
polygon(c(log10(xseq3), rev(log10(xseq3))), c(pred_ci_pcr_2[,3], rev(pred_ci_pcr_2[,4])), 
        col=adjustcolor(cols2[2], 0.3), border=NA)

lines(log10(xseq3), pred_ci_pcr_3[,1], lwd=3, col=cols2[1], lty=5)
polygon(c(log10(xseq3), rev(log10(xseq3))), c(pred_ci_pcr_3[,3], rev(pred_ci_pcr_3[,4])), 
        col=adjustcolor(cols2[1], 0.3), border=NA)

lines(log10(xseq3), pred_ci_pcr_4[,1], lwd=3, col=cols2[2], lty=5)
polygon(c(log10(xseq3), rev(log10(xseq3))), c(pred_ci_pcr_4[,3], rev(pred_ci_pcr_4[,4])), 
        col=adjustcolor(cols2[2], 0.3), border=NA)

axis(1, -1:8, c(10^(-1:3), format(10^(4:8), big.mark = ",", scientific = FALSE, trim=TRUE)),
     cex.axis = cx.axis)

legend("bottomright", c("HS-RDT, Uganda", "HS-RDT, Myanmar", "co-RDT, Uganda", "co-RDT, Myanmar"),
       col=rep(cols2, each = 2), lwd=3, lty=c(1,2,1,2), bty="n", cex = cx.leg)
mtext("D", side = 3, adj = -0.18, line=1, cex=1, font=2)

dev.off()


# calculate additional stats

bb1x = bb1[bb1$HRP2 >= 10 & bb1$HRP2 <= 1000, ]

table(bb1x$Pf.qRT.PCR..p.uL. > 0)
table(bb1x$Pf.qRT.PCR..p.uL. > 100)
table(bb1x$Pf.qRT.PCR..p.uL. > 1000)


bb2x = bb2[bb2$HRP2 >= 10 & bb2$HRP2 <= 1000, ]

table(bb2x$Parasitemia.qPCR > 0)
table(bb2x$Parasitemia.qPCR > 100)
table(bb2x$Parasitemia.qPCR > 1000)


# ALL CROSS-SECTIONAL DATA 
# PREVALENCE PLOTS
#install.packages("fmsb")
#install.packages("rmeta")

library(readxl)
library(dplyr)
library(binom)
library(rmeta)

dd = readxl::read_xlsx("../data from articles.xlsx", sheet = "cross-sectional", col_types = c("text", rep("numeric", 13)))

dd = data.frame(dd)

cols <- c("gray30", "gold2","deepskyblue", "firebrick2", "plum4", "darkorange1", "darkslategray3", "gray60","burlywood3", "yellow4", "deepskyblue3", "lightsalmon2", "mediumpurple4", "deeppink2", "yellow3", "brown4","mediumseagreen" , "slateblue2", "navy", "forestgreen")
#springgreen4

dd$col = cols


# plot 2 sensitivity of the HS-RDT

dd3 = dd[!is.na(dd$number.positiive.by.PCR),]
dd3 = arrange(dd3, desc(number.positiive.by.PCR))

dd3$pch_urdt = 16
dd3$pch_urdt[is.na(dd3$number.of.PCR.positives.that.are.also.RDT.positive)] = 17


cx.axis = 0.8
cx.names = 0.8
cx.leg = 0.6
cx.main = 0.8
cx = 0.6
cx.lab = 0.75

tiff("plots/FINAL_submitted/FIGURE_4.tiff", width = 170,height=104,units="mm",res=300, compression="lzw")

par(mar = c(4,4,2,12), xpd=FALSE)
plot(dd3$PCR.prevalence,
     dd3$number.of.PCR.positives.that.are.also.uRDT.positive/dd3$number.positiive.by.PCR*100,
     ylab = "Sensitivity (%)", xlab = "PCR prevalence (%)", las=1, ylim=c(0,100), col = dd3$col, pch=dd3$pch_urdt,lwd=2,
     cex=1.3, cex.lab = cx.lab, cex.axis = cx.axis)

mean_urdt = weighted.mean(dd3$number.of.PCR.positives.that.are.also.uRDT.positive/dd3$number.positiive.by.PCR*100,
                          dd3$number.positiive.by.PCR)

dd3_rdt = dd3[!is.na(dd3$number.positiive.by.RDT),]

mean_rdt = weighted.mean(dd3_rdt$number.of.PCR.positives.that.are.also.RDT.positive/dd3_rdt$number.positiive.by.PCR*100,
                         dd3_rdt$number.positiive.by.PCR)

# abline(h = mean_urdt, lwd=4, col="grey60")
# abline(h = c(mean_rdt+0.4, mean_rdt-0.4), lwd=1, col="grey60")
abline(h = seq(0, 100, by=10), lty=3, col="grey80")
abline(v = seq(0, 100, by=10), lty=3, col="grey80")

points(dd3$PCR.prevalence,
       dd3$number.of.PCR.positives.that.are.also.uRDT.positive/dd3$number.positiive.by.PCR*100,
       lwd=3, col=dd3$col, cex=1.3, pch=dd3$pch_urdt)

points(dd3$PCR.prevalence,
    dd3$number.of.PCR.positives.that.are.also.RDT.positive/dd3$number.positiive.by.PCR*100,
    lwd=3, col=dd3$col, cex=1.3)

segments(dd3$PCR.prevalence, 
         dd3$number.of.PCR.positives.that.are.also.RDT.positive/dd3$number.positiive.by.PCR*100,
         dd3$PCR.prevalence,
         dd3$number.of.PCR.positives.that.are.also.uRDT.positive/dd3$number.positiive.by.PCR*100)

dd3$urdt_sens = dd3$number.of.PCR.positives.that.are.also.RDT.positive/dd3$number.positiive.by.PCR*100

mm1 = glm(cbind(round(number.of.PCR.positives.that.are.also.uRDT.positive,0), 
          number.positiive.by.PCR - number.of.PCR.positives.that.are.also.uRDT.positive) ~ PCR.prevalence, data = dd3, family = "binomial")

summary(mm1)

xxseq = seq(0, 100, by=1)
preds = predict.glm(mm1, newdata = data.frame(PCR.prevalence = xxseq), se.fit = TRUE)

ilink <- family(mm1)$linkinv

lines(xxseq, ilink(preds$fit)*100, lwd=2, col="grey30")
polygon(c(xxseq, rev(xxseq)), 
        c(ilink(preds$fit + 2*preds$se.fit), rev(ilink(preds$fit - 2*preds$se.fit)))*100,
        border=NA, col = adjustcolor("grey", 0.4))


mm2 = glm(cbind(round(number.of.PCR.positives.that.are.also.RDT.positive,0), 
                number.positiive.by.PCR - number.of.PCR.positives.that.are.also.RDT.positive) ~ 
            PCR.prevalence, data = dd3, family = "binomial")

summary(mm2)

xxseq = seq(0, 100, by=1)
preds2 = predict.glm(mm2, newdata = data.frame(PCR.prevalence = xxseq), se.fit = TRUE)

ilink <- family(mm2)$linkinv

lines(xxseq, ilink(preds2$fit)*100, lwd=2, col="grey30", lty=5)
polygon(c(xxseq, rev(xxseq)), 
        c(ilink(preds2$fit + 2*preds2$se.fit), rev(ilink(preds2$fit - 2*preds2$se.fit)))*100,
        border=NA, col = adjustcolor("grey", 0.4))


dd3_fig = arrange(dd3, PCR.prevalence)

par(xpd=NA)
legend(86, 98, dd3_fig$Location, pch = dd3_fig$pch_urdt, col=dd3_fig$col, pt.cex = 1.1, title = "Study", bty="n",
       cex = cx)


legend(89, 4, c("HS-RDT estimated sensitivity", "Co-RDT estimated sensitivity"), col = "grey30", lwd=2,
       lty = c(1, 2), cex = cx)

# legend(89, 16, c(paste0("HS-RDT (",round(mean_urdt,1),"%)"), paste0("co-RDT (",round(mean_rdt,1),"%)")), col="grey60", lwd=c(4,5), title = "")
# legend(89, 16, c(paste0("HS-RDT (",round(mean_urdt,1),"%)"), paste0("co-RDT (",round(mean_rdt,1),"%)")), col=c(NA, "white"), lwd=c(4,3), bty="n", title = "Mean sensitivity")

dev.off()


# 
# weighted.ttest.ci <- function(x, weights, conf.level = 0.95) {
#   require(Hmisc)
#   nx <- length(x)
#   df <- nx - 1
#   vx <- wtd.var(x, weights, normwt = TRUE) ## From Hmisc
#   mx <- weighted.mean(x, weights)
#   stderr <- sqrt(vx/nx)
#   tstat <- mx/stderr ## not mx - mu
#   alpha <- 1 - conf.level
#   cint <- qt(1 - alpha/2, df)
#   cint <- tstat + c(-cint, cint)
#   cint * stderr
# }
# 
# 
# weighted.ttest.ci(x = dd3$number.of.PCR.positives.that.are.also.uRDT.positive/dd3$number.positiive.by.PCR*100,
#                   weights = dd3$number.positiive.by.PCR)
# 
# weighted.ttest.ci(x = dd3_rdt$number.of.PCR.positives.that.are.also.RDT.positive/dd3_rdt$number.positiive.by.PCR*100,
#                   weights = dd3_rdt$number.positiive.by.PCR)

library(readxl)
library(dplyr)
library(binom)
library(rmeta)

dd = readxl::read_xlsx("../data from articles.xlsx", sheet = "cross-sectional", col_types = c("text", rep("numeric", 13)))

dd = data.frame(dd)

cols <- c("gray30", "gold2","deepskyblue", "firebrick2", "plum4", "darkorange1", "darkslategray3", "gray60","burlywood3", "yellow4", "deepskyblue3", "lightsalmon2", "mediumpurple4", "deeppink2", "yellow3", "brown4","mediumseagreen" , "slateblue2",  "navy", "forestgreen")
#springgreen4

dd$col = cols

# calculate the confidence intervals
dd1 = dd[!is.na(dd$number.positiive.by.PCR),]

bc_pcr = binom.confint(dd1$number.positiive.by.PCR, dd1$number.tested.by.PCR, methods = "exact")
bc_pcr = data.frame(Study = dd1$Location, col = dd1$col, bc_pcr)

bc_urdt = binom.confint(dd1$number.positiive.by.uRDT, dd1$number.tested.by.uRDT, methods = "exact")
bc_urdt = data.frame(Study = dd1$Location, col = dd1$col, bc_urdt)

# Wu relationship
pcr_prev = seq(0, 1, l=101)

LO = function(x){log(x/(1-x))}
ULO = function(x){exp(x)/(1+exp(x))}
rdt_prev = ULO(-0.968 + 1.186*LO(pcr_prev))
sens = rdt_prev/pcr_prev

# fit this functional form to HS-RDT data

ff <- function(x){return(a0 + a1*LO(x))}

y = bc_urdt$mean
x = bc_pcr$mean

fit = nls(LO(y) ~ a0 + a1*LO(x), start = list(a0 = -1, a1 =1))

cx.axis = 0.8
cx.names = 0.8
cx.leg = 0.6
cx.main = 0.8
cx = 0.8
cx.lab = 0.75

tiff("plots/FINAL_submitted/FIGURE_3.tiff", width = 170,height=120,units="mm",res=300, compression="lzw")

par(mar=c(4,4,1.6,1), mfrow=c(1,2),oma=c(7.2,0,0,0), xpd=FALSE)
plot(bc_pcr$mean*100, bc_urdt$mean*100, pch=21, bg = "white", cex=cx,
     xlab = "PCR Prevalence (%)", ylab = "HS-RDT Prevalence (%)", las=1, xlim = c(0, 90), ylim=c(0,90),
     cex.axis = cx.axis, cex.lab = cx.lab)
abline(v = seq(0, 100, by=10), col="grey80", lty=3)
abline(h = seq(0, 100, by=10), col="grey80", lty=3)
rect(-2,-2,6,6, col=adjustcolor("pink", 0.7), border=NA)

abline(0,1, col="grey60", lty=5)
segments(bc_pcr$lower*100, bc_urdt$mean*100, bc_pcr$upper*100, bc_urdt$mean*100, col=as.character(bc_pcr$col), lwd=2)
segments(bc_pcr$mean*100, bc_urdt$lower*100, bc_pcr$mean*100, bc_urdt$upper*100, col=as.character(bc_pcr$col), lwd=2)
points(bc_pcr$mean*100, bc_urdt$mean*100, pch=21, bg = adjustcolor(bc_pcr$col, 0.8), cex=cx)

lines(pcr_prev*100, rdt_prev*100, col="orange", lwd=2, lty=5)

mtext("A", side = 3, adj=-0.1, font=2, cex=1.1)

par(mar=c(4,4,2.5,1))
plot(bc_pcr$mean*100, bc_urdt$mean*100, pch=21, bg = "white", cex=cx,
     xlab = "PCR Prevalence (%)", ylab = "HS-RDT Prevalence (%)", las=1, xlim = c(0, 6.5), ylim=c(0,6.5),
     cex.axis = cx.axis, cex.lab = cx.lab)
rect(-2,-2,6,6, col=adjustcolor("pink", 0.4), border=NA)
abline(0,1, col="grey60", lty=5)
abline(v = seq(0, 6, by=1), col="grey80", lty=3)
abline(h = seq(0, 6, by=1), col="grey80", lty=3)
segments(bc_pcr$lower*100, bc_urdt$mean*100, bc_pcr$upper*100, bc_urdt$mean*100, col=as.character(bc_pcr$col), lwd=2)
segments(bc_pcr$mean*100, bc_urdt$lower*100, bc_pcr$mean*100, bc_urdt$upper*100, col=as.character(bc_pcr$col), lwd=2)
points(bc_pcr$mean*100, bc_urdt$mean*100, pch=21, bg = adjustcolor(bc_pcr$col, 0.8), cex=cx)

lines(pcr_prev*100, rdt_prev*100, col="orange", lwd=2, lty=5)


mtext("B", side = 3, adj=-0.1, font=2, cex=1.1)

par(xpd=NA)

bc_pcr_leg = arrange(bc_pcr, mean)

legend(-12, -3.6, paste0(bc_pcr_leg$Study)[1:9], pch=21, 
       pt.bg = adjustcolor(bc_pcr_leg$col, 0.8)[1:9], pt.cex = 1.2, bty="n", cex=cx.leg)
legend(-4, -3.6, paste0(bc_pcr_leg$Study)[10:18], pch=21, 
       pt.bg = adjustcolor(bc_pcr_leg$col, 0.8)[10:18], pt.cex = 1.2, bty="n", cex=cx.leg)

legend(2, -3.6, "Fitted line from Wu\net al. PCR vs. co-RDT\nprevalence relationship", cex=cx.leg, col="orange", lwd=2, lty=2, bty="n")

dev.off()

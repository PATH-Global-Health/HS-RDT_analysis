rm(list = ls())

library(readxl)
library(dplyr)
library(binom)
library(rmeta)
library(fmsb)

dd = readxl::read_xlsx("../data from articles.xlsx", sheet = "cross-sectional", col_types = c("text", rep("numeric", 13)))

dd = data.frame(dd)

cols <- c("gray30", "gold2","deepskyblue", "firebrick2", "plum4", "darkorange1", "darkslategray3", "gray60","burlywood3", "yellow4", "deepskyblue3", "lightsalmon2", "mediumpurple4", "deeppink2", "yellow3", "brown4","mediumseagreen" , "slateblue2", "navy", "forestgreen")
#springgreen4

dd$col = cols

dd2 = dd[!is.na(dd$number.positiive.by.uRDT) & !is.na(dd$number.positiive.by.RDT),]


dd2$RDT.prevalence = dd2$number.positiive.by.RDT/dd2$number.tested.by.RDT
dd2$uRDT.prevalence = dd2$number.positiive.by.uRDT/dd2$number.tested.by.uRDT

dd2 = arrange(dd2, RDT.prevalence)


rr <- list()

for(i in 1:nrow(dd2)){
  rr[[i]] = riskratio(dd2$number.positiive.by.uRDT[i], dd2$number.positiive.by.RDT[i],
                      dd2$number.tested.by.uRDT[i], dd2$number.tested.by.RDT[i])
}

rr[[2]]$conf.int

getconfint <- function(x){ return(as.numeric(x$conf.int))}

rr2 = lapply(rr, getconfint)

cis = matrix(unlist(rr2), ncol=2, byrow = TRUE)

dd2$rr_lower = cis[,1]
dd2$rr_upper = cis[,2]

dd3=dd2
dd3 = dd2[-c(1),]

dd3 = arrange(dd3, dd3$uRDT.prevalence/dd3$RDT.prevalence)

dd3$Location[dd3$Location == "Haiti (Easy access groups exlcluding health facilties) (Druetz et al. 2020)"] = "Haiti (EAG exlcluding health facilties) (Druetz et al. 2020)"

tiff("plots/FINAL_submitted/FIGURE_5.tiff", width = 170,height=110,units="mm",res=300, compression="lzw")

par(mar = c(4,14,1,1))
plot(dd3$uRDT.prevalence/dd3$RDT.prevalence, nrow(dd3):1, xlim=c(0, 4), ylim = c(0, nrow(dd3)),
     xlab = "Ratio of HS-RDT prevalence to co-RDT prevalence",
     ylab="", yaxt="n", cex.lab = 0.7, cex.axis = 0.8)
abline(v = seq(0, 4, by=0.5), lty=3, col="grey70")
abline(v = 1, col="grey50", lty=5, lwd=2)
segments(dd3$rr_lower, nrow(dd3):1, dd3$rr_upper,nrow(dd3):1)
points(dd3$uRDT.prevalence/dd3$RDT.prevalence, nrow(dd3):1,pch=21,
       bg = dd3$col, cex=1.4)
axis(2, nrow(dd3):1, dd3$Location, las=1, cex.axis = 0.6)

meta = meta.DSL(number.tested.by.uRDT,number.tested.by.RDT,
               number.positiive.by.uRDT, number.positiive.by.RDT, data = dd3, statistic="RR")

# meta2 = meta.DSL(number.tested.by.uRDT,number.tested.by.RDT,
#                number.positiive.by.uRDT, number.positiive.by.RDT, data = dd3, statistic="RR")

polygon(c(1.26, 1.46, 1.7, 1.46), c(0, 0.2, 0, -0.2), col="deepskyblue1")

axis(2, 0, "Combined", las=1, cex.axis = 0.8)

dev.off()
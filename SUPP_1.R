
rm(list = ls())
library(readxl)
library(dplyr)
library(binom)
library(fmsb)
library(rmeta)
library(janitor)

dd = readxl::read_xlsx("../data from articles.xlsx", sheet = "cross-sectional", col_types = c("text", rep("numeric", 13)))

cols <- c("gray30", "gold2","deepskyblue", "firebrick2", "plum4", "darkorange1", "darkslategray3", "gray60","burlywood3", "yellow4", "deepskyblue3", "lightsalmon2", "mediumpurple4", "deeppink2", "yellow3", "brown4","mediumseagreen" , "slateblue2",  "navy", "forestgreen")

dd = dd %>% data.frame() %>% 
  clean_names() %>% 
  mutate(col = cols)



# PPV = number of true positives / (number of true positives + number of false positives)

dd$TP = dd$number_of_pcr_positives_that_are_also_u_rdt_positive 
dd$FP = dd$number_positiive_by_u_rdt - dd$number_of_pcr_positives_that_are_also_u_rdt_positive 

dd$PPV = dd$TP / (dd$TP + dd$FP)

dd$TP_co = dd$number_of_pcr_positives_that_are_also_rdt_positive 
dd$FP_co = dd$number_positiive_by_rdt - dd$number_of_pcr_positives_that_are_also_rdt_positive 

dd$PPV_co = dd$TP_co / (dd$TP_co + dd$FP_co)


#  NPV = number of true negatives / (number of true negatives + number of false negatives)



dd$FN = dd$number_positiive_by_pcr - dd$number_of_pcr_positives_that_are_also_u_rdt_positive

dd$TN = (dd$number_tested_by_u_rdt - dd$number_positiive_by_u_rdt) - dd$FN


dd$FN_co = dd$number_positiive_by_pcr - dd$number_of_pcr_positives_that_are_also_rdt_positive

dd$TN_co = (dd$number_tested_by_rdt - dd$number_positiive_by_rdt) - dd$FN_co

dd$NPV = dd$TN / (dd$TN + dd$FN)

dd$NPV_co = dd$TN_co / (dd$TN_co + dd$FN_co)

dd = dd %>% arrange(pcr_prevalence) 

bad = which(is.na(dd$PPV))

dd = dd[-bad, ]
dd = dd[-which(dd$location == "Uganda (Owalla et al. 2020)"),]

tiff("plots/FIG_SUPP_PPV.tiff", width = 240,height=140,units="mm",res=300, compression="lzw")

par(mar = c(4,4,2,20), xpd=FALSE)

plot(dd$pcr_prevalence, dd$PPV, ylim = c(0, 1), xlab = "PCR Prevalence (%)", ylab = "Positive predictive value", las=1,
     pch = 21, bg = adjustcolor(dd$col, 0.8), cex = 1.8)

points(dd$pcr_prevalence, dd$PPV_co, pch = 24, bg = adjustcolor(dd$col, 0.8))

segments(dd$pcr_prevalence, dd$PPV, dd$pcr_prevalence, dd$PPV_co)

par(xpd=NA)
legend(59, 1, dd$location, pch = 21, pt.bg  =adjustcolor(dd$col, 0.8), pt.cex = 1.4, title = "Study", bty="n", cex =1)

legend(59, 0, c("HS-RDT", "co-RDT"), pch = c(21, 24), pt.bg  = "white", pt.cex = 1.4, title = "RDT", bty="n", cex =1)

dev.off()


tiff("plots/FIG_SUPP_NPV.tiff", width = 240,height=140,units="mm",res=300, compression="lzw")

par(mar = c(4,4,2,20), xpd=FALSE)

plot(dd$pcr_prevalence, dd$NPV, ylim = c(0, 1), xlab = "PCR Prevalence (%)", ylab = "Negative predictive value", las=1,
     pch = 21, bg = adjustcolor(dd$col, 0.8), cex = 1.8)

points(dd$pcr_prevalence, dd$NPV_co, pch = 24, bg = adjustcolor(dd$col, 0.8))

segments(dd$pcr_prevalence, dd$NPV, dd$pcr_prevalence, dd$NPV_co)

par(xpd=NA)
legend(59, 1, dd$location, pch = 21, pt.bg  =adjustcolor(dd$col, 0.8), pt.cex = 1.4, title = "Study", bty="n", cex =1)

legend(59, 0, c("HS-RDT", "co-RDT"), pch = c(21, 24), pt.bg  = "white", pt.cex = 1.4, title = "RDT", bty="n", cex =1)

dev.off()


mean(dd$PPV)
weighted.mean(dd$PPV, dd$TP + dd$FP)
weighted.mean(dd$PPV, dd$TP)
weighted.mean(dd$PPV, dd$number_positiive_by_pcr)


mean(dd$PPV_co, na.rm=TRUE)
weighted.mean(dd$PPV_co, dd$number_positiive_by_pcr, na.rm = TRUE)

mean(dd$NPV)
weighted.mean(dd$NPV, dd$TN + dd$FN)
weighted.mean(dd$NPV, dd$TN)
weighted.mean(dd$NPV, dd$number_tested_by_pcr - dd$number_positiive_by_pcr)


mean(dd$NPV_co, na.rm=TRUE)
weighted.mean(dd$NPV_co, dd$number_tested_by_pcr - dd$number_positiive_by_pcr, na.rm=TRUE)

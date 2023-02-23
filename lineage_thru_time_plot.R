library("phytools")
library("phyloch")
library("strap")
library("coda")
library(TreeSim)
setwd("~/UAH_Docs/Fall_2021/Cave_beetles_proposal/Plot_geoscale_tree")
 
t <- read.beast("BEAST_DivTime_NoOutgroup_M0.0012SD0.059_combined_run1-2_summary_B10.tre")

plot(t,edge.color="red",cex=0.2)
ltt_plot <- ltt(t,log=FALSE)
lines(c(0, max(nodeHeights(t))), c(0, (length(t$tip.label))),lty = "dashed", lwd = 2, col = "red")
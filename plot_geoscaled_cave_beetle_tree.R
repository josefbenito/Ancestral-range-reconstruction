library("phytools")
library("phyloch")
library("strap")
library("coda")
setwd("~/UAH_Docs/Fall_2021/Cave_beetles_proposal/Plot_geoscale_tree")

t <- read.beast("BEAST_DivTime_NoOutgroup_M0.0012SD0.059_combined_run1-2_summary_B10.tre")
t$root.time <- t$height[1]

num_taxa <- length(t$tip.label)
display_all_node_bars <- FALSE

root_max <- t$"height_95%_HPD_MAX"[1]

pdffn = "Cave_beetle_DivTime_NoOutgroup_M0.0012SD0.059.pdf"
pdf(file=pdffn, width=8, height=9)

geoscalePhylo(tree=ladderize(t,right=FALSE), boxes="Age", cex.tip=1.0,cex.age=1.0,
				cex.ts=1.0,label.offset=0.2,x.lim=c(-15,root_max),lwd=1.5,width=1.5)

t$node.label<-t$posterior
p <- character(length(t$node.label))
p[t$node.label >= 0.95] <- "black"
p[t$node.label < 0.95 & t$node.label >= 0.75] <- "white"
p[t$node.label < 0.75] <- "red"
nodelabels(pch=21, cex=1.2, bg=p)

t$min_ages <- t$"height_95%_HPD_MIN"
t$max_ages <- t$"height_95%_HPD_MAX"

T1 <- get("last_plot.phylo", envir = .PlotPhyloEnv)

for(i in (Ntip(t) + 1):(t$Nnode + Ntip(t))) {
  lines(x = c(T1$root.time - t$min_ages[i - Ntip(t)],
              T1$root.time - t$max_ages[i - Ntip(t)]),
        y = rep(T1$yy[i], 2), lwd = 6, lend = 0,
        col = make.transparent("blue", 0.4))
}

dev.off()
cmdstr = paste0("open ", pdffn)
system(cmdstr)


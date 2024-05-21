# Format the pruned Bayesian time-calibrated phylogeny for BioGeoBEARS in tree format
library(ape)
library(phytools)
setwd('~/Google Drive/My Drive/projects/hirundinidae_phylogeny/workflow_reanalysis/analysis/biogeography/')

ha95 <- read.nexus('~/Google Drive/My Drive/projects/hirundinidae_phylogeny/workflow_reanalysis/analysis/beast/hirundinidae-all-rev1-95-prune.tree')
names <- read.table('hirundinidae-all-rev1-95-prune.short-taxon-names.txt',header=F)
ha95$tip.label <- names$V1
plot(ha95,cex=0.5)

write.tree(ha95,file='hirundinidae-all-rev1-95-prune.tree')


library(ape)
library(phytools)

outgroup <- c("Erythrocercus_mccallii_KU_8688", "Locustella_lanceolata_KU_4248")

### Concatenated maximum likelihood tree (ha95-rev1 dataset; n = 118 taxa)

# Read in and root tree; prune outgroups
ha95.rax <- read.tree("~/Google Drive/My Drive/projects/hirundinidae_phylogeny/workflow_reanalysis/analysis/concatenated/ha95-rev1/ha95-rev1-annotated.raxml.support.tree")
#write(ha95.rax$tip.label,'~/Google Drive/My Drive/projects/hirundinidae_phylogeny/workflow_reanalysis/plots/ha95-raxml-tips.txt') # Edit separate taxon and locality/clade labels from these below
ha95.rax <- root(ha95.rax, outgroup, edgelabel = T)
ha95.rax <- drop.tip(ha95.rax, outgroup)

# Format tip labels
taxon <- scan('~/Google Drive/My Drive/projects/hirundinidae_phylogeny/workflow_reanalysis/plots/ha95-raxml-tips-taxon.txt', what = "character")
locality <- scan('~/Google Drive/My Drive/projects/hirundinidae_phylogeny/workflow_reanalysis/plots/ha95-raxml-tips-locality.txt', what = "character")
taxon <- taxon[c(1:11, 14:118)] # Remove outgroup info; don't run these if outgroups are left in
locality <- locality[c(1:11, 14:118)]
ha95.rax$tip.label <- mixedFontLabel(taxon, locality, italic = 1)

# Format bootstrap support and organize by descending node order
boots.ha95.rax <- ha95.rax$node.label
boots.ha95.rax[boots.ha95.rax == "100"] <- ""
ha95.rax <- ladderize(ha95.rax)
plot(ha95.rax, cex = 0.7, no.margin = T, label.offset = 0.00025,edge.width = 1)
nodelabels(boots.ha95.rax, frame = "none", cex = 0.5, adj = -0.25)
add.scale.bar(cex = 0.7)

# Output at 14 x 10

### SVDquartets tree (ha95-rev1 dataset; n = 118 taxa)
ha95.svd <- read.nexus("~/Google Drive/My Drive/projects/hirundinidae_phylogeny/workflow_reanalysis/analysis/svdq/ha95-rev1/ha95-rev1.svdq.tre")
ha95.svd$edge.length <- NULL
ha95.svd <- root(ha95.svd, outgroup, edgelabel = T)
ha95.svd <- drop.tip(ha95.svd, outgroup)
boots.ha95.svd <- ha95.svd$node.label
boots.ha95.svd[boots.ha95.svd == "100"] <- ""
ha95.svd <- ladderize(ha95.svd)
plot(ha95.svd, cex = 0.6, no.margin = T, label.offset = 0.1, node.depth = 2)
nodelabels(boots.ha95.svd, frame= "none", cex = 0.5, adj = -0.25)

### RAxML tree (ha95-rev2 dataset; n = 103 taxa)
ha95.rax <- read.tree("~/Google Drive/My Drive/projects/hirundinidae_phylogeny/workflow_reanalysis/analysis/concatenated/ha95-rev2/ha95-rev2-annotated.raxml.support.tree")
ha95.rax <- root(ha95.rax, outgroup, edgelabel = T)
ha95.rax <- drop.tip(ha95.rax, outgroup)
boots.ha95.rax <- ha95.rax$node.label
boots.ha95.rax[boots.ha95.rax == "100"] <- ""
ha95.rax <- ladderize(ha95.rax)
plot(ha95.rax, cex = 0.7, no.margin = T, label.offset = 0.00025,edge.width = 1)
nodelabels(boots.ha95.rax, frame = "none", cex = 0.5, adj = -0.25)
add.scale.bar(cex = 0.7)

### SVDquartets tree (ha95-rev2 dataset; n = 103 taxa)
ha95.svd <- read.nexus("~/Google Drive/My Drive/projects/hirundinidae_phylogeny/workflow_reanalysis/analysis/svdq/ha95-rev2/ha95-rev2.svdq.tre")
ha95.svd$edge.length <- NULL
ha95.svd <- root(ha95.svd, outgroup, edgelabel = T)
ha95.svd <- drop.tip(ha95.svd, outgroup)
boots.ha95.svd <- ha95.svd$node.label
boots.ha95.svd[boots.ha95.svd == "100"] <- ""
ha95.svd <- ladderize(ha95.svd)
plot(ha95.svd, cex = 0.6, no.margin = T, label.offset = 0.1, node.depth = 2)
nodelabels(boots.ha95.svd, frame= "none", cex = 0.5, adj = -0.25)

### ASTRAL tree (ha95-rev3 dataset; n = 101 taxa)
ha95.astral <- read.tree("~/Google Drive/My Drive/projects/hirundinidae_phylogeny/workflow_reanalysis/analysis/astral/ha95-rev3/ha95-rev3.astral.tre")
ha95.astral$edge.length <- NULL
ha95.astral <- root(ha95.astral, outgroup, edgelabel = T)
ha95.astral <- drop.tip(ha95.astral, outgroup)
boots.ha95.astral <- ha95.astral$node.label
boots.ha95.astral[boots.ha95.astral == "1"] <- ""
ha95.astral <- ladderize(ha95.astral)
plot(ha95.astral, cex = 0.6, no.margin = T, label.offset = 0.1, node.depth = 2)
nodelabels(boots.ha95.astral, frame= "none", cex = 0.5, adj = -0.25)










#### OLDER STUFF BELOW----------------



### Compare RAxML and SVDquartets estimates (ha95-rev1 dataset; n = 118 taxa)
ha95.rax <- read.tree("~/Google Drive/My Drive/projects/hirundinidae_phylogeny/workflow_reanalysis/analysis/concatenated/ha95-rev1/ha95-rev1-annotated.raxml.support.tree")
ha95.rax <- root(ha95.rax, outgroup, edgelabel = T)
ha95.rax <- ladderize(ha95.rax)

ha95.svd <- read.nexus("~/Google Drive/My Drive/projects/hirundinidae_phylogeny/workflow_reanalysis/analysis/svdq/ha95-rev1/ha95-rev1.svdq.tre")
ha95.svd$edge.length <- NULL
ha95.svd <- root(ha95.svd, outgroup, edgelabel = T)
ha95.svd <- ladderize(ha95.svd)
plot(ha95.svd, cex = 0.6, no.margin = T, label.offset = 0.1, node.depth = 2)
nodelabels(ha95.svd$node.label, frame='none', cex = 0.5, adj = -0.25)


comparePhylo(ha95.rax, ha95.svd, plot = TRUE, force.rooted = TRUE,
             use.edge.length = FALSE, commons = FALSE,
             location = "bottomleft",cex=0.3)






### ASTRAL; all samples (skins + tissues), 95% complete matrix (Pilot Analysis)----
ral.t <- read.nexus("~/Google Drive/My Drive/projects/hirundinidae_phylogeny/workflow_reanalysis/results_pilot/species_tree/ha95-pilot/astral.ha95-pilot.50con.tre")
ral.t$edge.length <- NULL
ral.t <- root(ral.t, outgroup, edgelabel = T)
plot(ral.t, cex = 0.6, no.margin = T, label.offset = 0.1, node.depth = 2)

# ral.t <- drop.tip(ral.t, outgroup)

# Tips
write(ral.t$tip.label, "~/Desktop/tips.txt") # Edit them in Excel and text editor

genus.species <- scan("../../table_and_figure_components/astral-tissues-genus-species.txt",
                      what = "character")
locality <- scan("../../table_and_figure_components/astral-tissues-localities-etc.txt",
                 what = "character")

genus.species.no.out <- genus.species[c(1:62, 65:101)]
locality.no.out <- locality[c(1:62, 65:101)]

ral.t$tip.label <- mixedFontLabel(genus.species, locality, italic = 1)
# ral.t$tip.label <- mixedFontLabel(genus.species.no.out, locality.no.out, italic = 1)

# Boots - I should really fix the position of these. Maybe later.
boots.ral.t <- ral.t$node.label
boots.ral.t <- round(as.numeric(boots.ral.t)*100, 0)
boots.ral.t[boots.ral.t == "100"] <- ""
boots.ral.t[is.na(boots.ral.t)] <- ""

# Export
pdf(file = "~/Desktop/temp.pdf", height = 10, width = 8.5)
plot(ral.t, cex = 0.6, no.margin = T, label.offset = 0.1, node.depth = 2, edge.width = 1.25)
nodelabels(boots.ral.t, frame = "none", cex = 0.5, adj = -0.25)
dev.off()



### ORIGINAL CLARE CODE BELOW---------

#### RAxML all-V2 (annotated, skin branches truncated) [3.1] ####
rax.all.v2 <- read.tree("../raxml-ng/hirundinidae-all-95/ha95-annotated-truncated.raxml.support") 
rax.all.v2 <- root(rax.all.v2, outgroup, edgelabel = T)
rax.all.v2 <- drop.tip(rax.all.v2, outgroup)
plot(rax.all.v2, cex = 0.6, no.margin = T, label.offset = 0.0005)
add.scale.bar()

# Add a root edge
rax.all.v2$root.edge <- 0.001
plot(rax.all.v2, cex = 0.6, no.margin = T, label.offset = 0.0001, root.edge = T)

# Tips
write(rax.all.v2$tip.label, "~/Desktop/tips.txt") # Edit them in Excel and text editor

genus.species <- scan("../../table_and_figure_components/raxml-all-v2-genus-species.txt",
                      what = "character")
locality <- scan("../../table_and_figure_components/raxml-all-v2-localites-etc.txt", what = "character")
rax.all.v2$tip.label <- mixedFontLabel(genus.species, locality, italic = 1)

# Boots - I should really fix the position of these. Maybe later.
boots.rax.all.v2 <- rax.all.v2$node.label
boots.rax.all.v2[boots.rax.all.v2 == "100"] <- ""

# Export
pdf(file = "~/Desktop/temp.pdf", height = 10, width = 8)
plot(rax.all.v2, cex = 0.6, no.margin = T, root.edge = T, edge.width = 1.25, label.offset = 0.0001)
add.scale.bar(x = 0.002, y = 0, cex = 0.8)
nodelabels(boots.rax.all.v2, frame = "none", cex = 0.5, adj = -0.25)
dev.off()


#### ASTRAL, tissues only [3.2] ####
ral.t <- read.nexus("../speciestrees/astral-ht95.50con.tre")
ral.t$edge.length <- NULL
ral.t <- root(ral.t, outgroup, edgelabel = T)
plot(ral.t, cex = 0.6, no.margin = T, label.offset = 0.1, node.depth = 2)

# ral.t <- drop.tip(ral.t, outgroup)

# Tips
write(ral.t$tip.label, "~/Desktop/tips.txt") # Edit them in Excel and text editor

genus.species <- scan("../../table_and_figure_components/astral-tissues-genus-species.txt",
                      what = "character")
locality <- scan("../../table_and_figure_components/astral-tissues-localities-etc.txt",
                 what = "character")

genus.species.no.out <- genus.species[c(1:62, 65:101)]
locality.no.out <- locality[c(1:62, 65:101)]

ral.t$tip.label <- mixedFontLabel(genus.species, locality, italic = 1)
# ral.t$tip.label <- mixedFontLabel(genus.species.no.out, locality.no.out, italic = 1)

# Boots - I should really fix the position of these. Maybe later.
boots.ral.t <- ral.t$node.label
boots.ral.t <- round(as.numeric(boots.ral.t)*100, 0)
boots.ral.t[boots.ral.t == "100"] <- ""
boots.ral.t[is.na(boots.ral.t)] <- ""

# Export
pdf(file = "~/Desktop/temp.pdf", height = 10, width = 8.5)
plot(ral.t, cex = 0.6, no.margin = T, label.offset = 0.1, node.depth = 2, edge.width = 1.25)
nodelabels(boots.ral.t, frame = "none", cex = 0.5, adj = -0.25)
dev.off()


#### Supplementary trees ####

#### RAxML tissues (annotated) [S3.2] ####
rax.t <- read.tree("../raxml-ng/hirundinidae-tissues-95/ht95-annotated.raxml.support") 
rax.t <- root(rax.t, outgroup, edgelabel = T)
rax.t <- drop.tip(rax.t, outgroup)
plot(rax.t, cex = 0.6, no.margin = T, label.offset = 0.0002)
add.scale.bar()

# Add a root edge
rax.t$root.edge <- 0.001
plot(rax.t, cex = 0.6, no.margin = T, label.offset = 0.0002, root.edge = T)

# Boots - I should really fix the position of these. Maybe later.
boots.t <- rax.t$node.label
boots.t[boots.t == "100"] <- ""

# Export
pdf(file = "~/Desktop/temp.pdf", height = 10, width = 8)
plot(rax.t, cex = 0.6, no.margin = T, root.edge = T, edge.width = 1.25, label.offset = 0.0002)
add.scale.bar(x = 0.002, y = 0, cex = 0.8)
nodelabels(boots.t, frame = "none", cex = 0.5, adj = -0.25)
dev.off()


#### ASTRAL, all-v3 only [3.3] ####
ral.all.v3 <- read.nexus("../speciestrees/astral-ha-v3-95.50con.tre")
ral.all.v3$edge.length <- NULL
ral.all.v3 <- root(ral.all.v3, outgroup, edgelabel = T)
plot(ral.all.v3, cex = 0.6, no.margin = T, label.offset = 0.1, node.depth = 2)

# Boots - I should really fix the position of these. Maybe later.
boots.all.v3 <- ral.all.v3$node.label
boots.all.v3 <- round(as.numeric(boots.all.v3)*100, 0)
boots.all.v3[boots.all.v3 == "100"] <- ""
boots.all.v3[is.na(boots.all.v3)] <- ""

# Export
pdf(file = "~/Desktop/temp.pdf", height = 10, width = 8.5)
plot(ral.all.v3, cex = 0.6, no.margin = T, label.offset = 0.1, node.depth = 2, edge.width = 1.25)
nodelabels(boots.all.v3, frame = "none", cex = 0.5, adj = -0.25)
dev.off()


#### ASTRID, tissues only [3.4] ####
rid.t <- read.nexus("../speciestrees/astrid-ht95.50con.tre")
rid.t$edge.length <- NULL
rid.t <- root(rid.t, outgroup, edgelabel = T)
plot(rid.t, cex = 0.6, no.margin = T, label.offset = 0.1, node.depth = 2)

# Boots - I should really fix the position of these. Maybe later.
boots.rid.t <- rid.t$node.label
boots.rid.t <- round(as.numeric(boots.rid.t)*100, 0)
boots.rid.t[boots.rid.t == "100"] <- ""
boots.rid.t[is.na(boots.rid.t)] <- ""

# Export
pdf(file = "~/Desktop/temp.pdf", height = 10, width = 8.5)
plot(rid.t, cex = 0.6, no.margin = T, label.offset = 0.1, node.depth = 2, edge.width = 1.25)
nodelabels(boots.rid.t, frame = "none", cex = 0.5, adj = -0.25)
dev.off()


#### ASTRID, all-v3 [3.5] ####
rid.all.v3 <- read.nexus("../speciestrees/astrid-ha-v3-95.50con.tre")
rid.all.v3$edge.length <- NULL
rid.all.v3 <- root(rid.all.v3, outgroup, edgelabel = T)
plot(rid.all.v3, cex = 0.6, no.margin = T, label.offset = 0.1, node.depth = 2)

# Boots - I should really fix the position of these. Maybe later.
boots.rid.all.v3 <- rid.all.v3$node.label
boots.rid.all.v3 <- round(as.numeric(boots.rid.all.v3)*100, 0)
boots.rid.all.v3[boots.rid.all.v3 == "100"] <- ""
boots.rid.all.v3[is.na(boots.rid.all.v3)] <- ""

# Export
pdf(file = "~/Desktop/temp.pdf", height = 10, width = 8.5)
plot(rid.all.v3, cex = 0.6, no.margin = T, label.offset = 0.1, node.depth = 2, edge.width = 1.25)
nodelabels(boots.rid.all.v3, frame = "none", cex = 0.5, adj = -0.25)
dev.off()




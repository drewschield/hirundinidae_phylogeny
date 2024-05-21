library(ape)
library(phytools)
library(diversitree)
library(geiger)

# Set the working directory
wd <- "~/Google Drive/My Drive/projects/hirundinidae_phylogeny/workflow_reanalysis/analysis/nest_sociality_reconstructions/"
setwd(wd)
getwd()

#### Load the swallow tree ####
# This is the RAxML tree inferred from the 95 percent complete data matrix of all samples,
# pruned to 81 taxa and transformed into a time tree.
tr <- read.tree("./hirundinidae-all-rev1-95-prune.tree")
plot(tr, cex = 0.6, no.margin = F, label.offset = 0.1); axisPhylo()

length(tr$tip.label) # 88
tr$Nnode #87

# Drop Pseudochelidon sirintarae
tr <- drop.tip(tr, "Pseudochelidon_sirintarae")
plot(tr, cex = 0.6, no.margin = F, label.offset = 0.1); axisPhylo()


#### Load and prepare the nest type data ####
# Nest type dataset for the 80 taxa in the tree, tip names should match.
# Nest type as a six-state character.
  # excavate burrow, adopt burrow, adopt cavity, mud open cup, mud closed cup, mud retort
dat1 <- read.delim("./swallow_data_nest.txt", stringsAsFactors = F)
names(dat1)
head(dat1)
  # "nest" is the verbal coding
  # "taxon" is self explanatory
nest <- dat1[,1:2]
summary(as.factor(nest$nest))


# Merge the burrow adopters and cavity adopters
nest$nest[nest$nest == "adopt_burrow"] <- "adopt_cavity"
summary(as.factor(nest$nest))

# Get the data frame and the tree into the same order, to make plotting easier later.
dim(nest); length(tr$tip.label); all(tr$tip.label %in% nest$taxon) # The tips and taxon names match...
all.equal(tr$tip.label, nest$taxon) # ...but are not in the same order 
rownames(nest) <- nest$taxon
nest <- nest[tr$tip.label,] # get the data frame in the same order as the tree
all.equal(tr$tip.label, nest$taxon) # Now the tree and the data frame are in the same order.

# convert the nest data to a binary matrix
nest.types <- c("adopt_cavity","excavate_burrow","mud_enclosed_cup","mud_open_cup","mud_retort")
nest.vec <- nest$nest
names(nest.vec) <- rownames(nest)
nest.mat <- to.matrix(nest.vec, nest.types)


#### Load and prepare the sociality data ####
# Breeding sociality data for the 80 taxa in the tree, tip names should match.
# Breeding sociality as a three-state character.
# solitary, solitary and small group, colony
  # six taxa are uncertain, and are coded as solitary/small_group
dat2 <- read.delim("./swallow_data_sociality.txt", stringsAsFactors = F)
names(dat2)
head(dat2)
  # "breeding_sociality" is the sociality category
  # "taxon" is self explanatory
soc <- dat2
summary(as.factor(soc$breeding_sociality))

# Get the data frame and the tree into the same order, to make plotting easier later.
dim(soc); length(tr$tip.label); all(tr$tip.label %in% soc$taxon) # The tips and taxon names match...
all.equal(tr$tip.label, soc$taxon) # ...but are not in the same order 
rownames(soc) <- soc$taxon
soc <- soc[tr$tip.label,] # get the data frame in the same order as the tree
all.equal(tr$tip.label, soc$taxon) # Now the tree and the data frame are in the same order.

# Convert the sociality data to a binary matrix
soc.types <- c("colony","small_group","solitary")
soc.vec <- soc$breeding_sociality
names(soc.vec) <- rownames(soc)
soc.mat <- to.matrix(soc.vec, soc.types)

# Fix the uncertain taxa - give them 50/50 probabilities.
temp <- (soc$taxon[soc$breeding_sociality == "solitary/small_group"])
for (i in 1:length(temp)) {
  soc.mat[temp[i],] <- c(0,0.5,0.5)
}


#### Plot nest type and sociality onto the tree ####

# Nest type alone
nest.cols <- c("#abd9e9","#2c7bb6","#fdae61","#ffffbf","#d7191c")
plot(tr, no.margin = T, cex = 0.7, label.offset = 0.5)
tiplabels(pie = nest.mat, piecol = nest.cols, cex = 0.3)

# Sociality alone
soc.cols <- c("#636363","#bdbdbd","#f0f0f0") 
plot(tr, no.margin = T, cex = 0.7, label.offset = 0.5)
tiplabels(pie = soc.mat, piecol = soc.cols, cex = 0.3)

# Both characters together
pdf("./temp.pdf", height = 11, width = 8.5)
par(mar = c(0,2,0,0), xpd=T)
plot(tr, no.margin = F, cex = 0.8, label.offset = 1.3, edge.width = 1.5)
tiplabels(pie = nest.mat, piecol = nest.cols, cex = 0.35)
tiplabels(pie = soc.mat, piecol = soc.cols, cex = 0.35, offset = 0.75)
legend(x = -2, y = 0, xjust = 0, yjust = 0,
       title = "BREEDING SOCIALITY", title.adj = 0,
       legend = c("solitary","solitary and small groups","colonial"),
       pch = 21, pt.bg = c("#f0f0f0","#bdbdbd","#636363"), pt.cex = 1.5, bty = "n", cex = 0.8)
legend(x = -2, y = 7, xjust = 0, yjust = 0,
       title = "NEST TYPE", title.adj = 0,
       legend = c("excavates burrow","adopts burrow or cavity",
                  "mud nest: open cup","mud nest: enclosed cup","mud nest: retort"),
       pch = 21, pt.bg = c("#2c7bb6","#abd9e9","#ffffbf","#fdae61","#d7191c"),
       pt.cex = 1.5, bty = "n", cex = 0.8)
dev.off()


#### Nest type evolution: model comparisons - fitDiscrete{geiger} ####

#### 1. Simplified data matrix: which transitions can be eliminated from the model? ####

# Reduce mud nest construction to one state
nest.3st.vec <- nest.vec
nest.3st.vec[nest.3st.vec == "mud_open_cup"] <- "mud_nest"
nest.3st.vec[nest.3st.vec == "mud_enclosed_cup"] <- "mud_nest"
nest.3st.vec[nest.3st.vec == "mud_retort"] <- "mud_nest"

summary(as.factor(nest.3st.vec))

# Fit and compare models of nest type evolution, nest type as three-state character.
results.3st.fitDiscrete <- data.frame(model=c("Model1.ER","Model2.SYM","Model3.ARD","Model4","Model5"),
                                      k=numeric(5), lnL=numeric(5), AICc=numeric(5),
                                      deltaAICc=numeric(5), Akaike_weight=numeric(5))

model1 <- fitDiscrete(tr, nest.3st.vec, model = "ER")
results.3st.fitDiscrete[1,-1] <- c(model1$opt$k, model1$opt$lnL, model1$opt$aicc, 0, 0)

model2 <- fitDiscrete(tr, nest.3st.vec, model = "SYM")
results.3st.fitDiscrete[2,-1] <- c(model2$opt$k, model2$opt$lnL, model2$opt$aicc, 0, 0)

model3 <- fitDiscrete(tr, nest.3st.vec, model = "ARD")
results.3st.fitDiscrete[3,-1]<- c(model3$opt$k, model3$opt$lnL, model3$opt$aicc, 0, 0)

model4 <- fitDiscrete(tr, nest.3st.vec, model = matrix(c(NA,1,2,0,NA,0,3,4,NA),3,3))
results.3st.fitDiscrete[4,-1]<- c(model4$opt$k, model4$opt$lnL, model4$opt$aicc, 0, 0)

model5 <- fitDiscrete(tr, nest.3st.vec, model = matrix(c(NA,1,2,0,NA,0,0,3,NA),3,3))
results.3st.fitDiscrete[5,-1]<- c(c(model5$opt$k, model5$opt$lnL, model5$opt$aicc, 0, 0))

# # Save models 1-5
save(model1, file = "fitDiscrete_model1.RData")
save(model2, file = "fitDiscrete_model2.RData")
save(model3, file = "fitDiscrete_model3.RData")
save(model4, file = "fitDiscrete_model4.RData")
save(model5, file = "fitDiscrete_model5.RData")

# load("fitDiscrete_model1.RData")
# load("fitDiscrete_model2.RData")
# load("fitDiscrete_model3.RData")
# load("fitDiscrete_model4.RData")
# load("fitDiscrete_model5.RData")

# deltaAICc and Akaike weight
results.3st.fitDiscrete$deltaAICc <- results.3st.fitDiscrete$AICc - min(results.3st.fitDiscrete$AICc)
results.3st.fitDiscrete$Akaike_weight <-
  exp(-0.5 * results.3st.fitDiscrete$deltaAICc) / sum(exp(-0.5 * results.3st.fitDiscrete$deltaAICc))

results.3st.fitDiscrete <- results.3st.fitDiscrete[order(results.3st.fitDiscrete$AICc), ]
results.3st.fitDiscrete

write.table(cbind(results.3st.fitDiscrete[ , 1:2], round(results.3st.fitDiscrete[ , 3:6], 3)),
            file = "fitDiscrete_model_comparison_3state.txt", quote = F, sep = "\t")

# I could also do a likelihood ratio test of Model4 vs. Model5.
# Model4 (3 rates) is nested within Model5 (4 rates):
  # Model5 adds a parameter to Model4: cavity adoption --> mud nest construction
# LRT H0: the addition of the fourth rate to Model5 does not increase the likelihood of the data
# LRT Ha: the addition of the fourth rate to Model5 increases the likelihood of the data
# Given that these two models have identical likelihoods, I actually already know that I can't
# reject the null - the difference between the two likelihoods is zero.

# I'll compare the ARD (model3) and SYM (model2) models instead.
teststat <- 2*(model3$opt$lnL - model2$opt$lnL)
pchisq(teststat, df = (model3$opt$k - model2$opt$k), lower.tail = F)
  # [1]  0.078
  # So, I can kind of marginally reject the null. The additional parameters (rates) in the ARD
  # model sort of confer a greater likelihood on the nest data.


#### 2. Fitting 5-state models: Which type of mud-nest arose first? ####

# Models 6-8: 4 parameters
  # all mud --> adopt cavity equal, all rates among mud nests equal
model6 <-
  fitDiscrete(tr, nest.vec,
              model = rbind(c(NA,0,0,0,0), c(1,NA,0,4,0), c(2,0,NA,3,3), c(2,0,3,NA,3), c(2,0,3,3,NA)))
model7 <-
  fitDiscrete(tr, nest.vec,
              model = rbind(c(NA,0,0,0,0), c(1,NA,3,0,0), c(2,0,NA,4,4), c(2,0,4,NA,4), c(2,0,4,4,NA)))
model8 <-
  fitDiscrete(tr, nest.vec,
              model = rbind(c(NA,0,0,0,0), c(1,NA,0,0,4), c(2,0,NA,3,3), c(2,0,3,NA,3), c(2,0,3,3,NA)))

# Model 9: 6 parameters
  # Does not specify which mud nest type arose first
model9 <-
  fitDiscrete(tr, nest.vec,
              model = rbind(c(NA,0,0,0,0), c(1,NA,3,5,6), c(2,0,NA,4,4), c(2,0,4,NA,4), c(2,0,4,4,NA)))


# # Save models 6-8
save(model6, file = "fitDiscrete_model6.RData")
save(model7, file = "fitDiscrete_model7.RData")
save(model8, file = "fitDiscrete_model8.RData")
save(model9, file = "fitDiscrete_model9.RData")

# load("fitDiscrete_model6.RData")
# load("fitDiscrete_model7.RData")
# load("fitDiscrete_model8.RData")
# load("fitDiscrete_model9.RData")

# Compare models 6-8
models.5st <- list(model6,model7,model8)

results.5st.fitDiscrete <- data.frame(model=6:8, hypothesis=c("open","enclosed","retort"),
                                      k=numeric(3), lnL=numeric(3), AICc=numeric(3),
                                      deltaAICc=numeric(3) ,Akaike_weight=numeric(3))

for (i in 1:length(models.5st)) {
  results.5st.fitDiscrete[i, c(-1,-2)] <-
    c(models.5st[[i]]$opt$k, models.5st[[i]]$opt$lnL, models.5st[[i]]$opt$aicc, 0, 0) 
}

# deltaAICc and Akaike weight
results.5st.fitDiscrete$deltaAICc <- results.5st.fitDiscrete$AICc - min(results.5st.fitDiscrete$AICc)
results.5st.fitDiscrete$Akaike_weight <-
  exp(-0.5 * results.5st.fitDiscrete$deltaAICc) / sum(exp(-0.5 * results.5st.fitDiscrete$deltaAICc))

results.5st.fitDiscrete <- results.5st.fitDiscrete[order(results.5st.fitDiscrete$AICc), ]
results.5st.fitDiscrete

write.table(cbind(results.5st.fitDiscrete[ , 1:3], round(results.5st.fitDiscrete[ , 4:7], 3)),
            file = "fitDiscrete_model_comparison_5state.txt",
            quote=FALSE, sep="\t")

# Look at model 9 too.
model9


#### 3. Ancestral character estimation: Diversitree ####
## Model 6 (open cup first) is fairly strongly preferred, but I'll run ACE with all four models
## just to see how they look.
## I am calculating marginal ancestral state probabilities here: asr.marginal()

## Prepare the data.
nest.types # defined above
nest.vec.code <- nest.vec
for (i in 1:length(nest.types)) {
  nest.vec.code[nest.vec.code == nest.types[i]] <- i
}
nest.vec.code <- as.numeric(nest.vec.code); names(nest.vec.code) <- names(nest.vec)


## Make the likelihood function for this tree and dataset.
  # mkn = multistate Markov model

lik.mkn <- make.mkn(tr, nest.vec.code, k=5, control = list(root=ROOT.FLAT))
argnames(lik.mkn)

## Model 6: open mud cup first - best model according to AICc, Akaike weight

# constrain the model
lik.mkn.m6 <- constrain(lik.mkn,
                        q12 ~ 0, q13 ~ 0, q14 ~ 0, q15 ~ 0, q23 ~ 0, q25 ~ 0, q32 ~ 0, q42 ~ 0, q52 ~ 0,
                        q41 ~ q31, q51 ~ q31,
                        q35 ~ q34, q43 ~ q34, q45 ~ q34, q53 ~ q34, q54 ~ q34)
argnames(lik.mkn.m6)

# get a starting point for the ML search
p.mkn <- starting.point.musse(tr, 5)

# fit the model
fitted.mkn.m6 <- find.mle(lik.mkn.m6, p.mkn[argnames(lik.mkn.m6)])

# model parameters - the four estimated transition rates, and the log likelihood
fitted.mkn.m6[1:2]

# get the marginal ancestral states for the nodes
  # returns a matrix of ancestral state probabilities for the nodes
  # Note that: asr.marginal(lik.mkn.m6, coef(fitted.mkn.m6)) returns the ancestral state probabilities
  # with the states rows and the nodes as columns, transposing makes it much easier to plot (and see).
anc.states.m6 <- t(asr.marginal(lik.mkn.m6, coef(fitted.mkn.m6)))

# plot the ancestral states
par(mar = c(0,2,0,0), xpd=T)
plot(tr, no.margin = F, cex = 0.8, label.offset = 0.5, edge.width = 1.5)
tiplabels(pie = nest.mat, piecol = nest.cols, cex = 0.35)
nodelabels(pie = anc.states.m6, piecol = nest.cols, cex = 0.4)
arrows(3.25,62,4.3,59.5, col = "red", lwd = 2, angle = 20, length = 0.15)
legend(x = -2.25, y = 0, xjust = 0, yjust = 0,
       title = "NEST TYPE", title.adj = 0,
       legend = c("excavates burrow","adopts burrow or cavity",
                  "mud nest: open cup","mud nest: enclosed cup","mud nest: retort"),
       pch = 21, pt.bg = c("#2c7bb6","#abd9e9","#ffffbf","#fdae61","#d7191c"),
       pt.cex = 1.5, bty = "n", cex = 0.8)

#### Plot model6 (best model) to pdf, with sociality ####

# plot the ancestral states, include sociality as a second tip state
pdf("./temp.pdf", height = 11, width = 8.5)
par(mar = c(0,2,0,0), xpd=T)
plot(tr, no.margin = F, cex = 0.8, label.offset = 1.3, edge.width = 1.5)
tiplabels(pie = nest.mat, piecol = nest.cols, cex = 0.35)
tiplabels(pie = soc.mat, piecol = soc.cols, cex = 0.35, offset = 0.75)
nodelabels(pie = anc.states.m6, piecol = nest.cols, cex = 0.4)
arrows(3.25,62,4.3,59.5, col = "red", lwd = 2, angle = 20, length = 0.15)
legend(x = -2.25, y = 0, xjust = 0, yjust = 0,
       title = "BREEDING SOCIALITY", title.adj = 0,
       legend = c("solitary","solitary and small groups","colonial"),
       pch = 21, pt.bg = c("#f0f0f0","#bdbdbd","#636363"), pt.cex = 1.5, bty = "n", cex = 0.8)
legend(x = -2.25, y = 7, xjust = 0, yjust = 0,
       title = "NEST TYPE", title.adj = 0,
       legend = c("excavates burrow","adopts burrow or cavity",
                  "mud nest: open cup","mud nest: enclosed cup","mud nest: retort"),
       pch = 21, pt.bg = c("#2c7bb6","#abd9e9","#ffffbf","#fdae61","#d7191c"),
       pt.cex = 1.5, bty = "n", cex = 0.8)
dev.off()

## Model 7: enclosed mud cup first

# constrain the model
lik.mkn.m7 <- constrain(lik.mkn,
                        q12 ~ 0, q13 ~ 0, q14 ~ 0, q15 ~ 0, q24 ~ 0, q25 ~ 0, q32 ~ 0, q42 ~ 0, q52 ~ 0,
                        q41 ~ q31, q51 ~ q31,
                        q35 ~ q34, q43 ~ q34, q45 ~ q34, q53 ~ q34, q54 ~ q34)
argnames(lik.mkn.m7)

# get a starting point for the ML search
p.mkn <- starting.point.musse(tr, 5)

# fit the model
fitted.mkn.m7 <- find.mle(lik.mkn.m7, p.mkn[argnames(lik.mkn.m7)])

# model parameters - the four estimated transition rates, and the log likelihood
fitted.mkn.m7[1:2]

# get the marginal ancestral states for the nodes
# returns a matrix of ancestral state probabilities for the nodes
anc.states.m7 <- t(asr.marginal(lik.mkn.m7, coef(fitted.mkn.m7)))

# plot the ancestral states
par(mar = c(0,2,0,0), xpd=T)
plot(tr, no.margin = F, cex = 0.8, label.offset = 0.5, edge.width = 1.5)
tiplabels(pie = nest.mat, piecol = nest.cols, cex = 0.35)
nodelabels(pie = anc.states.m7, piecol = nest.cols, cex = 0.4)
arrows(3.25,62,4.3,59.5, col = "red", lwd = 2, angle = 20, length = 0.15)
legend(x = -2.25, y = 0, xjust = 0, yjust = 0,
       title = "NEST TYPE", title.adj = 0,
       legend = c("excavates burrow","adopts burrow or cavity",
                  "mud nest: open cup","mud nest: enclosed cup","mud nest: retort"),
       pch = 21, pt.bg = c("#2c7bb6","#abd9e9","#ffffbf","#fdae61","#d7191c"),
       pt.cex = 1.5, bty = "n", cex = 0.8)


## Model 8: mud retort first

# constrain the model
lik.mkn.m8 <- constrain(lik.mkn,
                       q12 ~ 0, q13 ~ 0, q14 ~ 0, q15 ~ 0, q23 ~ 0, q24 ~ 0, q32 ~ 0, q42 ~ 0, q52 ~ 0,
                       q41 ~ q31, q51 ~ q31,
                       q35 ~ q34, q43 ~ q34, q45 ~ q34, q53 ~ q34, q54 ~ q34)
argnames(lik.mkn.m8)

# get a starting point for the ML search
p.mkn <- starting.point.musse(tr, 5)

# fit the model
fitted.mkn.m8 <- find.mle(lik.mkn.m8, p.mkn[argnames(lik.mkn.m8)])

# model parameters - the four estimated transition rates, and the log likelihood
fitted.mkn.m8[1:2]

# get the marginal ancestral states for the nodes
# returns a matrix of ancestral state probabilities for the nodes
anc.states.m8 <- t(asr.marginal(lik.mkn.m8, coef(fitted.mkn.m8)))

# plot the ancestral states
par(mar = c(0,2,0,0), xpd=T)
plot(tr, no.margin = F, cex = 0.8, label.offset = 0.5, edge.width = 1.5)
tiplabels(pie = nest.mat, piecol = nest.cols, cex = 0.35)
nodelabels(pie = anc.states.m8, piecol = nest.cols, cex = 0.4)
arrows(3.25,62,4.3,59.5, col = "red", lwd = 2, angle = 20, length = 0.15)
legend(x = -2.25, y = 0, xjust = 0, yjust = 0,
       title = "NEST TYPE", title.adj = 0,
       legend = c("excavates burrow","adopts burrow or cavity",
                  "mud nest: open cup","mud nest: enclosed cup","mud nest: retort"),
       pch = 21, pt.bg = c("#2c7bb6","#abd9e9","#ffffbf","#fdae61","#d7191c"),
       pt.cex = 1.5, bty = "n", cex = 0.8)


## Model 9: no a priori assumption about which mud nest arose first

# constrain the model
lik.mkn.m9 <- constrain(lik.mkn,
                        q12 ~ 0, q13 ~ 0, q14 ~ 0, q15 ~ 0, q32 ~ 0, q42 ~ 0, q52 ~ 0,
                        q41 ~ q31, q51 ~ q31,
                        q35 ~ q34, q43 ~ q34, q45 ~ q34, q53 ~ q34, q54 ~ q34)
argnames(lik.mkn.m9)

# get a starting point for the ML search
p.mkn <- starting.point.musse(tr, 5)

# fit the model
fitted.mkn.m9 <- find.mle(lik.mkn.m9, p.mkn[argnames(lik.mkn.m9)])

# model parameters - the four estimated transition rates, and the log likelihood
fitted.mkn.m9[1:2]

# get the marginal ancestral states for the nodes
# returns a matrix of ancestral state probabilities for the nodes
anc.states.m9 <- t(asr.marginal(lik.mkn.m9, coef(fitted.mkn.m9)))

# plot the ancestral states
par(mar = c(0,2,0,0), xpd=T)
plot(tr, no.margin = F, cex = 0.8, label.offset = 0.5, edge.width = 1.5)
tiplabels(pie = nest.mat, piecol = nest.cols, cex = 0.35)
nodelabels(pie = anc.states.m9, piecol = nest.cols, cex = 0.4)
arrows(3.25,62,4.3,59.5, col = "red", lwd = 2, angle = 20, length = 0.15)
legend(x = -2.25, y = 0, xjust = 0, yjust = 0,
       title = "NEST TYPE", title.adj = 0,
       legend = c("excavates burrow","adopts burrow or cavity",
                  "mud nest: open cup","mud nest: enclosed cup","mud nest: retort"),
       pch = 21, pt.bg = c("#2c7bb6","#abd9e9","#ffffbf","#fdae61","#d7191c"),
       pt.cex = 1.5, bty = "n", cex = 0.8)


# Models 6 and 9 return almost exactly identical ancestral state probabilities:

rbind(100*round(anc.states.m6[41,],4), 100*round(anc.states.m9[41,],4)) # ancestor of all mud nesters
rbind(100*round(anc.states.m6[42,],4), 100*round(anc.states.m9[42,],4)) # ancestor of enclosed/retort


# Plot Models 6 and 9 together. 
# See: http://blog.phytools.org/2018/02/a-refresher-on-plotting-facing.html

#tips.original <- tr$tip.label 
#tr$tip.label <- c("A._pileata","A._tibialis","A._fasciata","P._cyanoleuca","P._melanoleuca","A_fucata",
#                  "O._murina","O._andecola","O._flavipes","S._serripennis_ridgwayi","S._serripennis",
#                  "S._ruficollis","P._tapera","P._subis","P._elegans","P._dominicensis","P._sinaloae",
#                  "P._chalybea","P._murphyi","T._bicolor","T._cyaneoviridis","T._euchrysea",
#                  "T._thalassina_brachyptera","T._thalassina_thalassina","T._leucorrhoa",
#                  "T._leucopyga","T._stolzmanni","T._albiventer","T._albilinea","R._riparia.",
#                  "R._diluta","R._paludicola_cowani","R._paludicola_paludicola","P._borbonica",
#                  "P._brazzae","N._cincta","C._leucosterna","P._griseopyga","P._fluvicola","P._ariel",
#                  "P._nigricans","P._pyrrhonota","P._fulva","P._rufocollaris","P._preussi",
#                  "P._spilodera","C._abyssinica","C._cucullata","C._daurica_Eur._SE_Asia",
#                  "C._daurica_Africa","C._senegalensis","C._semirufa_semirufa","C._semirufa_gordoni",
#                  "P._fuliginosa","D._nipalense","D._dasypus","D._urbicum_lagopodum",
#                  "D._urbicum_urbicum","H._atrocaerulea","H._dimidiata","H._rustica","H._angolensis",
#                  "H._lucida","H._aethiopica","H._albigularis","H._nigrita","H._smithii","H._neoxena",
#                  "H._tahitica","P._fuligula_S_pop.","P._fuligula_N_pop.","P._concolor","P._rupestris",
#                  "P._pristoptera_N","P._albiceps","P._pristoptera_S","P._fuliginosa","P._obscura",
#                  "P._nitens","P._eurystomina")

pdf("./temp.pdf", height = 8.5, width = 11)

layout(matrix(1:3,1,3),widths=c(0.45,0.1,0.45))
plotTree(tr,ftype="off",ylim=c(-1,Ntip(tr)))
nodelabels(pie = anc.states.m6, piecol = nest.cols, cex = 0.7)
tiplabels(pie = nest.mat, piecol = nest.cols, cex = 0.45)
title(main = "MODEL 6", line = -1.5)
legend(x = 0, y = 0, xjust = 0, yjust = 0,
      title = "NEST TYPE", title.adj = 0,
      legend = c("excavates burrow","adopts burrow or cavity",
                 "mud nest: open cup","mud nest: enclosed cup","mud nest: retort"),
      pch = 21, pt.bg = c("#2c7bb6","#abd9e9","#ffffbf","#fdae61","#d7191c"),
      pt.cex = 2.5, bty = "n", cex = 1.25)


ylim<-get("last_plot.phylo",envir=.PlotPhyloEnv)$y.lim
plot.new(); plot.window(xlim=c(-0.1,0.1),ylim=ylim)
text(rep(0,Ntip(tr)),1:Ntip(tr),gsub("_"," ",tr$tip.label),cex=0.9,
     font=3)

plotTree(tr,ftype="off",ylim=c(-1,Ntip(tr)),direction="leftwards")
nodelabels(pie = anc.states.m9, piecol = nest.cols, cex = 0.7)
tiplabels(pie = nest.mat, piecol = nest.cols, cex = 0.45)
title(main = "MODEL 9", line = -1.5)

dev.off()

tr$tip.label <- tips.original


#### 4. Ancestral character estimation: stochastic character mapping (phytools) ####

# Mixed model of nest type evolution: COMBINED STOCHASTIC CHARACTER MAPPING
# simulate a number of character histories proportionate to Akaike weights
# Derive Akaike weights and transition rates from simmap - fit models

# Functions to calculate AIC and Akaike weight
aic <- function(logL, k) 2*k-2*logL
aic.w<-function(aic) {
  d.aic <- aic-min(aic)
  exp(-1/2*d.aic)/sum(exp(-1/2*d.aic))
}

# Make a new object for the data, which should be a named vector of character states.
x <- nest.vec

# Models of nest type evolution
nest.types # matrix order should match order of this vector
m6 <- rbind(c(0,0,0,0,0), c(1,0,0,4,0), c(2,0,0,3,3), c(2,0,3,0,3), c(2,0,3,3,0)) # open
m7 <- rbind(c(0,0,0,0,0), c(1,0,3,0,0), c(2,0,0,4,4), c(2,0,4,0,4), c(2,0,4,4,0)) # enclosed
m8 <- rbind(c(0,0,0,0,0), c(1,0,0,0,4), c(2,0,0,3,3), c(2,0,3,0,3), c(2,0,3,3,0)) # retort

models <- list(m6,m7,m8) # put them in a list

# Number of model parameters
k <- rep(4,3)


## Get the likelihoods for the models, fitted to the swallow tree and nest data with make.simmap

L <- rep(0, length(models))
names(L) <- c("m6","m7","m8")
for (i in 1:length(models)) {
  L[i] <- make.simmap(tr, x, model = models[[i]], pi="equal", Q="empirical")$logL 
}
L # save as a separate object


## Calculate AIC, Akaike weight, and nsim for each model
# nsim = number of character histories to simulate under each model of nest-type evolution

AIC <- mapply(aic, L, k); AIC
AICw <- aic.w(AIC); AICw
AICw <- round(AICw,3); AICw

round(cbind(L,AIC,AICw), 3)

nsim <- 1000 # Simulate 1000 stochastic character histories in total

# Calculate how many character histories to simulate under each model of nest-type evolution
Nsim <- nsim*round(AICw,3) # round AICw to get an integer when multiplied by 1000
d <- if(sum(Nsim) > nsim) -1 else 1
nsim <- Nsim + d*sample(c(rep(1, abs(nsim-sum(Nsim))),
                      rep(0, length(Nsim) -abs(nsim-sum(Nsim)))))

nsim<-nsim[nsim!=0] #remove any with nsim==0

round(cbind(L,AIC,AICw,nsim), 3)

## Simulate the character histories with make.simmap, Q = "empirical"

# This code simulates trees for all models at once, outputting them as a single list. 
# for (i in 1:length(nsim)) {
#   obj <- make.simmap(tr, x, model=models[[i]], nsim=nsim[i],
#                      pi="equal", Q="empirical")
#   if (nsim[i]==1) {
#     obj <- list(obj)
#     class(obj) <- "multiPhylo"
#   }
#   if (i==1) {
#     trees <- obj
#   } else {
#     trees <- c(trees, obj)
#   }
# }
# save(trees, file = "scm_combined_trees_empirical.RData")

# It would be nice to be able to look at the trees for the different models independently.
# This code creates a separate object for the character histories simulated under each character
# evolution model, and then combines then.

## Simulate the character histories with make.simmap
# Q = "empirical"

# Open mud nest first: m6
if (nsim["m6"]!=0) {
  scm.m6 <- make.simmap(tr, x, model=m6, nsim=nsim["m6"], pi="equal", Q="empirical")
  if(nsim["m6"]==1) {
    scm.m6 <- list(scm.m6)
    class(scm.m6) <- "multiPhylo"
  }
}
# Enclosed mud nest first: m7
if (nsim["m7"]!=0) {
  scm.m7 <- make.simmap(tr, x, model=m7, nsim=nsim["m7"], pi="equal", Q="empirical")
  if(nsim["m7"]==1) {
    scm.m7 <- list(scm.m7)
    class(scm.m7) <- "multiPhylo"
  }
}
# Mud retort first: m8
if (nsim["m8"]!=0) {
  scm.m8 <- make.simmap(tr, x, model=m8, nsim=nsim["m8"], pi="equal", Q="empirical")
  if(nsim["m8"]==1) {
    scm.m8 <- list(scm.m8)
    class(scm.m8) <- "multiPhylo"
  }
}

save(scm.m6, file = "scm_trees/scm_empirical_m6_trees.RData")
save(scm.m7, file = "scm_trees/scm_empirical_m7_trees.RData")
save(scm.m8, file = "scm_trees/scm_empirical_m8_trees.RData")

# load("scm_trees/scm_empirical_m6_trees.RData")
# load("scm_trees/scm_empirical_m7_trees.RData")
# load("scm_trees/scm_empirical_m8_trees.RData")

## Simulate the character histories with make.simmap
# Q = "mcmc"

# Open mud nest first: m6
if (nsim["m6"]!=0) {
  scm.m6 <- make.simmap(tr, x, model=m6, nsim=nsim["m6"], pi="equal", Q="mcmc")
  if(nsim["m6"]==1) {
    scm.m6 <- list(scm.m6)
    class(scm.m6) <- "multiPhylo"
  }
}
# Enclosed mud nest first: m7
if (nsim["m7"]!=0) {
  scm.m7 <- make.simmap(tr, x, model=m7, nsim=nsim["m7"], pi="equal", Q="mcmc")
  if(nsim["m7"]==1) {
    scm.m7 <- list(scm.m7)
    class(scm.m7) <- "multiPhylo"
  }
}
# Mud retort first: m8
if (nsim["m8"]!=0) {
  scm.m8 <- make.simmap(tr, x, model=m8, nsim=nsim["m8"], pi="equal", Q="mcmc")
  if(nsim["m8"]==1) {
    scm.m8 <- list(scm.m8)
    class(scm.m8) <- "multiPhylo"
  }
}

# save(scm.m6, file = "scm_trees/scm_mcmc_m6_trees.RData")
# save(scm.m7, file = "scm_trees/scm_mcmc_m7_trees.RData")
# save(scm.m8, file = "scm_trees/scm_mcmc_m8_trees.RData")

# load("scm_trees/scm_mcmc_m6_trees.RData")
# load("scm_trees/scm_mcmc_m7_trees.RData")
# load("scm_trees/scm_mcmc_m8_trees.RData")


# Combine the mapped histories
scm.combined <- c(scm.m6, scm.m7, scm.m8)


## Look at the results

# Summarize the posterior density
pd.combined <- describe.simmap(scm.combined)
pd.m6 <- describe.simmap(scm.m6); pd.m7 <- describe.simmap(scm.m7); pd.m8 <- describe.simmap(scm.m8)

# Get the number of transitions among character states (nest types)
trans <- rbind(mean=colMeans(pd.combined$count),
               min=apply(pd.combined$count,2,min), max=apply(pd.combined$count,2,max))

# I'm interested in transitions from burrow excavation to each of the mud nest types
  # I only expect one transition per tree
colnames(trans) # columns 7-9, the first column summarizes all 21 types of transitions

trans <- trans[,c(1,8,7,9)]
colnames(trans) <- c("N","excavate --> open","excavate --> enclosed","excavate --> retort")
  # Note that the mean number of transitions exactly recapitulates the original AIC weights.
trans


# Plot with pie charts
pdf("~/Desktop/temp_Qmcmc.pdf", height = 11, width = 8.5)
par(mar = c(0,2,0,0), xpd=T)
plot(tr, no.margin = F, cex = 0.8, label.offset = 0.7, edge.width = 1.5)
tiplabels(pie = nest.mat, piecol = nest.cols, cex = 0.35)
legend(x = -2.25, y = 0, xjust = 0, yjust = 0,
       title = "NEST TYPE", title.adj = 0,
       legend = c("excavates burrow","adopts burrow or cavity",
                  "mud nest: open cup","mud nest: enclosed cup","mud nest: retort"),
       pch = 21, pt.bg = c("#2c7bb6","#abd9e9","#ffffbf","#fdae61","#d7191c"),
       pt.cex = 1.5, bty = "n", cex = 0.8)
nodelabels(node=1:tr$Nnode+Ntip(tr), pie=pd.combined$ace, piecol=nest.cols, cex=0.4)
dev.off()

# I can't do a density map, because my character has more than two states.


#### 5. Analyses with Petrochelidon fuliginosa categorized as building a retort-type mud nest ####
## Model testing and ancestral character estimation with Diversitree

# Change Petrochelidon fuliginosa nest type to "mud_retort"
nest.vec["Petrochelidon_fuliginosa"] <- "mud_retort"


## Fit the 5-state models: Which type of mud-nest arose first?

# Models 6-8: 4 parameters
# all mud --> adopt cavity equal, all rates among mud nests equal
model6 <-
  fitDiscrete(tr, nest.vec,
              model = rbind(c(NA,0,0,0,0), c(1,NA,0,4,0), c(2,0,NA,3,3), c(2,0,3,NA,3), c(2,0,3,3,NA)))
model7 <-
  fitDiscrete(tr, nest.vec,
              model = rbind(c(NA,0,0,0,0), c(1,NA,3,0,0), c(2,0,NA,4,4), c(2,0,4,NA,4), c(2,0,4,4,NA)))
model8 <-
  fitDiscrete(tr, nest.vec,
              model = rbind(c(NA,0,0,0,0), c(1,NA,0,0,4), c(2,0,NA,3,3), c(2,0,3,NA,3), c(2,0,3,3,NA)))

# Model 9: 6 parameters
# Does not specify which mud nest type arose first
model9 <-
  fitDiscrete(tr, nest.vec,
              model = rbind(c(NA,0,0,0,0), c(1,NA,3,5,6), c(2,0,NA,4,4), c(2,0,4,NA,4), c(2,0,4,4,NA)))


# Compare models 6-8
models.5st <- list(model6,model7,model8)

results.5st.fitDiscrete <- data.frame(model=6:8, hypothesis=c("open","enclosed","retort"),
                                      k=numeric(3), lnL=numeric(3), AICc=numeric(3),
                                      deltaAICc=numeric(3) ,Akaike_weight=numeric(3))

for (i in 1:length(models.5st)) {
  results.5st.fitDiscrete[i, c(-1,-2)] <-
    c(models.5st[[i]]$opt$k, models.5st[[i]]$opt$lnL, models.5st[[i]]$opt$aicc, 0, 0) 
}

# deltaAICc and Akaike weight
results.5st.fitDiscrete$deltaAICc <- results.5st.fitDiscrete$AICc - min(results.5st.fitDiscrete$AICc)
results.5st.fitDiscrete$Akaike_weight <-
  exp(-0.5 * results.5st.fitDiscrete$deltaAICc) / sum(exp(-0.5 * results.5st.fitDiscrete$deltaAICc))

results.5st.fitDiscrete <- results.5st.fitDiscrete[order(results.5st.fitDiscrete$AICc), ]
results.5st.fitDiscrete
cbind(results.5st.fitDiscrete[ , 1:3], round(results.5st.fitDiscrete[ , 4:7], 3))

# Look at model 9 too.
model9

## Ancestral character estimation with Diversitree, under the unconstrained model (Model 9)
# I am calculating marginal ancestral state probabilities here: asr.marginal()

## Prepare the data.
nest.types # defined above
nest.vec.code <- nest.vec
for (i in 1:length(nest.types)) {
  nest.vec.code[nest.vec.code == nest.types[i]] <- i
}
nest.vec.code <- as.numeric(nest.vec.code); names(nest.vec.code) <- names(nest.vec)


## Make the likelihood function for this tree and dataset.
# mkn = multistate Markov model

lik.mkn <- make.mkn(tr, nest.vec.code, k=5, control = list(root=ROOT.FLAT))
argnames(lik.mkn)

## Model 9: no a priori assumption about which mud nest arose first

# constrain the model
lik.mkn.m9 <- constrain(lik.mkn,
                        q12 ~ 0, q13 ~ 0, q14 ~ 0, q15 ~ 0, q32 ~ 0, q42 ~ 0, q52 ~ 0,
                        q41 ~ q31, q51 ~ q31,
                        q35 ~ q34, q43 ~ q34, q45 ~ q34, q53 ~ q34, q54 ~ q34)
argnames(lik.mkn.m9)

# get a starting point for the ML search
p.mkn <- starting.point.musse(tr, 5)

# fit the model
fitted.mkn.m9 <- find.mle(lik.mkn.m9, p.mkn[argnames(lik.mkn.m9)])

# model parameters - the four estimated transition rates, and the log likelihood
fitted.mkn.m9[1:2]

# get the marginal ancestral states for the nodes
# returns a matrix of ancestral state probabilities for the nodes
anc.states.m9 <- t(asr.marginal(lik.mkn.m9, coef(fitted.mkn.m9)))

# plot the ancestral states
par(mar = c(0,2,0,0), xpd=T)
plot(tr, no.margin = F, cex = 0.8, label.offset = 0.7, edge.width = 1.5)
tiplabels(pie = nest.mat, piecol = nest.cols, cex = 0.35)
nodelabels(pie = anc.states.m9, piecol = nest.cols, cex = 0.4)
legend(x = -2.25, y = 0, xjust = 0, yjust = 0,
       title = "NEST TYPE", title.adj = 0,
       legend = c("excavates burrow","adopts burrow or cavity",
                  "mud nest: open cup","mud nest: enclosed cup","mud nest: retort"),
       pch = 21, pt.bg = c("#2c7bb6","#abd9e9","#ffffbf","#fdae61","#d7191c"),
       pt.cex = 1.5, bty = "n", cex = 0.8)


# Ancestral state probabilities
rbind(nest.types, 100*round(anc.states.m9[41,],4)) # ancestor of all mud nesters
rbind(nest.types, 100*round(anc.states.m9[42,],4)) # ancestor of retor and enclosed nest


#### 6. Analyses with more Delichon tips, Delichon ancestor pushed back ####
## Model testing and ancestral character estimation with Diversitree

# In order to keep the tree ultrametric, all numbers have to have the same number of digits.
options(digits=10)

# This is the clade of Delichon + Petrochelidon fuliginosa from the ultrametric tree.
tr <- read.tree(text = "(Petrochelidon_fuliginosa:6.532304003,((Delichon_nipalense:2.948680035,(Delichon_dasypus:2.538025704,Delichon_urbicum_lagopodum:2.538025704)100:0.4106543312)100:0.3897331553,Delichon_urbicum_urbicum:3.338413191)100:3.193890812)100:2.232046653;")
plot(tr)

# Halve the length of the branch to Delichon
3.193890812/2 # 1.596945406
tr <- read.tree(text = "(Petrochelidon_fuliginosa:6.532304003,((Delichon_nipalense:2.948680035,(Delichon_dasypus:2.538025704,Delichon_urbicum_lagopodum:2.538025704)100:0.4106543312)100:0.3897331553,Delichon_urbicum_urbicum:3.338413191)100:1.596945406)100:2.232046653;")
plot(tr)

# And add that length to two of the Delichon edges
tr$edge.length # edges 3 to 8
for (i in c(3,8)) {
  tr$edge.length[i] <- tr$edge.length[i] + 1.596945406
}
plot(tr)
is.ultrametric(tr)

# Now split the Delichon tips
tr$tip.label <- c("Petrochelidon_fuliginosa","A","B","C","D")
plot(tr)
print(write.tree(tr))

# A
2.948680035-0.75; # 2.198680035
tr <- read.tree(text = "(Petrochelidon_fuliginosa:6.532304003,(((Aa:2.198680035,Ab:2.198680035):0.75,(B:2.538025704,C:2.538025704)100:0.4106543312)100:1.986678561,D:4.935358597)100:1.596945406)100:2.232046653;")
plot(tr)

# B
2.538025704-0.75 # 1.788025704
tr <- read.tree(text = "(Petrochelidon_fuliginosa:6.532304003,(((Aa:2.198680035,Ab:2.198680035):0.75,((Ba:1.788025704,Bb:1.788025704):0.75,C:2.538025704)100:0.4106543312)100:1.986678561,D:4.935358597)100:1.596945406)100:2.232046653;")
plot(tr)

# C
2.538025704-1.0 # 1.538025704
tr <- read.tree(text = "(Petrochelidon_fuliginosa:6.532304003,(((Aa:2.198680035,Ab:2.198680035):0.75,((Ba:1.788025704,Bb:1.788025704):0.75,(Ca:1.538025704,Cb:1.538025704):1.0)100:0.4106543312)100:1.986678561,D:4.935358597)100:1.596945406)100:2.232046653;")
plot(tr)

# D
4.935358597-0.75 # 4.185358597
tr <- read.tree(text = "(Petrochelidon_fuliginosa:6.532304003,(((Aa:2.198680035,Ab:2.198680035):0.75,((Ba:1.788025704,Bb:1.788025704):0.75,(Ca:1.538025704,Cb:1.538025704):1.0)100:0.4106543312)100:1.986678561,(Da:4.185358597,Db:4.185358597):0.75)100:1.596945406)100:2.232046653;")
plot(tr)

# Da, Db
4.185358597-0.75 # 3.435358597
4.185358597-1.0 # 3.185358597
tr <- read.tree(text = "(Petrochelidon_fuliginosa:6.532304003,(((Aa:2.198680035,Ab:2.198680035):0.75,((Ba:1.788025704,Bb:1.788025704):0.75,(Ca:1.538025704,Cb:1.538025704):1.0)100:0.4106543312)100:1.986678561,((Daa:3.435358597,Dab:3.435358597):0.75,(Dba:3.185358597,Dbb:3.185358597):1.0):0.75)100:1.596945406)100:2.232046653;")
plot(tr)

# Aa
2.198680035-0.5 # 1.698680035
tr <- read.tree(text = "(Petrochelidon_fuliginosa:6.532304003,((((Aaa:1.698680035,Aab:1.698680035):0.5,Ab:2.198680035):0.75,((Ba:1.788025704,Bb:1.788025704):0.75,(Ca:1.538025704,Cb:1.538025704):1.0)100:0.4106543312)100:1.986678561,((Daa:3.435358597,Dab:3.435358597):0.75,(Dba:3.185358597,Dbb:3.185358597):1.0):0.75)100:1.596945406)100:2.232046653;")
plot(tr)

# Daa
3.435358597-0.75 # 2.685358597
tr <- read.tree(text = "(Petrochelidon_fuliginosa:6.532304003,((((Aaa:1.698680035,Aab:1.698680035):0.5,Ab:2.198680035):0.75,((Ba:1.788025704,Bb:1.788025704):0.75,(Ca:1.538025704,Cb:1.538025704):1.0)100:0.4106543312)100:1.986678561,(((Daaa:2.685358597,Daab:2.685358597):0.75,Dab:3.435358597):0.75,(Dba:3.185358597,Dbb:3.185358597):1.0):0.75)100:1.596945406)100:2.232046653;")
plot(tr)

# Rename tips
tr$tip.label <- c("Petrochelidon_fuliginosa", paste(rep("Delichon",12), LETTERS[1:12], sep = "_"))
plot(tr)

is.ultrametric(tr)

print(write.tree(tr)) # Paste this text into the chronogram newick file.

### And proceed with analyses:

## Load the modified tree
# Delichon now has 12 tips
tr <- read.tree("chr-ha95-pruned-modified-delichon.tre")
plot(tr, cex = 0.6, no.margin = T, label.offset = 0.1)

# Drop Pseudochelidon sirintarae
tr <- drop.tip(tr, "Pseudochelidon_sirintarae")

## Load, modify, and prepare the nest type data
dat1 <- read.delim("swallow_nest_data_tree_taxa.txt", stringsAsFactors = F)
names(dat1)
head(dat1)
nest <- dat1[,1:2]
summary(as.factor(nest$nest))

# Merge the burrow adopters and cavity adopters
nest$nest[nest$nest == "adopt_burrow"] <- "adopt_cavity"
summary(as.factor(nest$nest))

# Remove the real Delichon species and add the artificial taxa, all as enclosed nest builders.
tr$tip.label # artificial Delichon are tips 55-66
delichon <- data.frame(taxon=tr$tip.label[55:66], nest=rep("mud_enclosed_cup", 12), stringsAsFactors = F)

nest$taxon # rows 13-16
nest <- nest[c(1:12,17:80),]
nest <- rbind(nest, delichon)


# Get the data frame and the tree into the same order, to make plotting easier later.
dim(nest); length(tr$tip.label); all(tr$tip.label %in% nest$taxon) # The tips and taxon names match...
all.equal(tr$tip.label, nest$taxon) # ...but are not in the same order 
rownames(nest) <- nest$taxon
nest <- nest[tr$tip.label,] # get the data frame in the same order as the tree
all.equal(tr$tip.label, nest$taxon) # Now the tree and the data frame are in the same order.

# The mud nesters are now much more evenly divided among types
summary(as.factor(nest$nest))

# convert the nest data to a binary matrix
nest.types <- c("adopt_cavity","excavate_burrow","mud_enclosed_cup","mud_open_cup","mud_retort")
nest.vec <- nest$nest
names(nest.vec) <- rownames(nest)
nest.mat <- to.matrix(nest.vec, nest.types)


## Plot nest type onto the tree

# Nest type alone
nest.cols <- c("#abd9e9","#2c7bb6","#fdae61","#ffffbf","#d7191c")
plot(tr, no.margin = T, cex = 0.7, label.offset = 0.5)
tiplabels(pie = nest.mat, piecol = nest.cols, cex = 0.3)


## Fit the 5-state models: Which type of mud-nest arose first?

# Models 6-8: 4 parameters
# all mud --> adopt cavity equal, all rates among mud nests equal
model6 <-
  fitDiscrete(tr, nest.vec,
              model = rbind(c(NA,0,0,0,0), c(1,NA,0,4,0), c(2,0,NA,3,3), c(2,0,3,NA,3), c(2,0,3,3,NA)))
model7 <-
  fitDiscrete(tr, nest.vec,
              model = rbind(c(NA,0,0,0,0), c(1,NA,3,0,0), c(2,0,NA,4,4), c(2,0,4,NA,4), c(2,0,4,4,NA)))
model8 <-
  fitDiscrete(tr, nest.vec,
              model = rbind(c(NA,0,0,0,0), c(1,NA,0,0,4), c(2,0,NA,3,3), c(2,0,3,NA,3), c(2,0,3,3,NA)))

# Model 9: 6 parameters
# Does not specify which mud nest type arose first
model9 <-
  fitDiscrete(tr, nest.vec,
              model = rbind(c(NA,0,0,0,0), c(1,NA,3,5,6), c(2,0,NA,4,4), c(2,0,4,NA,4), c(2,0,4,4,NA)))


# Compare models 6-8
models.5st <- list(model6,model7,model8)

results.5st.fitDiscrete <- data.frame(model=6:8,
                                      hypothesis=c("open","enclosed","retort"),
                                      k=numeric(3), lnL=numeric(3), AICc=numeric(3),
                                      deltaAICc=numeric(3) ,Akaike_weight=numeric(3))

for (i in 1:length(models.5st)) {
  results.5st.fitDiscrete[i, c(-1,-2)] <-
    c(models.5st[[i]]$opt$k, models.5st[[i]]$opt$lnL, models.5st[[i]]$opt$aicc, 0, 0) 
}

# deltaAICc and Akaike weight
results.5st.fitDiscrete$deltaAICc <- results.5st.fitDiscrete$AICc - min(results.5st.fitDiscrete$AICc)
results.5st.fitDiscrete$Akaike_weight <-
  exp(-0.5 * results.5st.fitDiscrete$deltaAICc) / sum(exp(-0.5 * results.5st.fitDiscrete$deltaAICc))

results.5st.fitDiscrete <- results.5st.fitDiscrete[order(results.5st.fitDiscrete$AICc), ]
results.5st.fitDiscrete
cbind(results.5st.fitDiscrete[ , 1:3], round(results.5st.fitDiscrete[ , 4:7], 3))

# Look at model 9 too.
model9

## Ancestral character estimation with Diversitree, Model 6

## Prepare the data.
nest.types # defined above
nest.vec.code <- nest.vec
for (i in 1:length(nest.types)) {
  nest.vec.code[nest.vec.code == nest.types[i]] <- i
}
nest.vec.code <- as.numeric(nest.vec.code); names(nest.vec.code) <- names(nest.vec)


## Make the likelihood function for this tree and dataset.
# mkn = multistate Markov model

lik.mkn <- make.mkn(tr, nest.vec.code, k=5, control = list(root=ROOT.FLAT))
argnames(lik.mkn)

## Model 6: open mud cup first

# constrain the model
lik.mkn.m6 <- constrain(lik.mkn,
                        q12 ~ 0, q13 ~ 0, q14 ~ 0, q15 ~ 0, q23 ~ 0, q25 ~ 0, q32 ~ 0, q42 ~ 0, q52 ~ 0,
                        q41 ~ q31, q51 ~ q31,
                        q35 ~ q34, q43 ~ q34, q45 ~ q34, q53 ~ q34, q54 ~ q34)
argnames(lik.mkn.m6)

# get a starting point for the ML search
p.mkn <- starting.point.musse(tr, 5)

# fit the model
fitted.mkn.m6 <- find.mle(lik.mkn.m6, p.mkn[argnames(lik.mkn.m6)])

# model parameters - the four estimated transition rates, and the log likelihood
fitted.mkn.m6[1:2]

# get the marginal ancestral states for the nodes
# returns a matrix of ancestral state probabilities for the nodes
anc.states.m6 <- t(asr.marginal(lik.mkn.m6, coef(fitted.mkn.m6)))

# plot the ancestral states
par(mar = c(0,2,0,0), xpd=T)
plot(tr, no.margin = F, cex = 0.8, label.offset = 0.5, edge.width = 1.5)
tiplabels(pie = nest.mat, piecol = nest.cols, cex = 0.35)
nodelabels(pie = anc.states.m6, piecol = nest.cols, cex = 0.4)
legend(x = -3.5, y = 0, xjust = 0, yjust = 0,
       title = "NEST TYPE", title.adj = 0,
       legend = c("excavates burrow","adopts burrow or cavity",
                  "mud nest: open cup","mud nest: enclosed cup","mud nest: retort"),
       pch = 21, pt.bg = c("#2c7bb6","#abd9e9","#ffffbf","#fdae61","#d7191c"),
       pt.cex = 1.5, bty = "n", cex = 0.8)


## Model 9: no a priori assumption about which mud nest arose first

# constrain the model
lik.mkn.m9 <- constrain(lik.mkn,
                        q12 ~ 0, q13 ~ 0, q14 ~ 0, q15 ~ 0, q32 ~ 0, q42 ~ 0, q52 ~ 0,
                        q41 ~ q31, q51 ~ q31,
                        q35 ~ q34, q43 ~ q34, q45 ~ q34, q53 ~ q34, q54 ~ q34)
argnames(lik.mkn.m9)

# get a starting point for the ML search
p.mkn <- starting.point.musse(tr, 5)

# fit the model
fitted.mkn.m9 <- find.mle(lik.mkn.m9, p.mkn[argnames(lik.mkn.m9)])

# model parameters - the four estimated transition rates, and the log likelihood
fitted.mkn.m9[1:2]

# get the marginal ancestral states for the nodes
# returns a matrix of ancestral state probabilities for the nodes
anc.states.m9 <- t(asr.marginal(lik.mkn.m9, coef(fitted.mkn.m9)))

# plot the ancestral states
par(mar = c(0,2,0,0), xpd=T)
plot(tr, no.margin = F, cex = 0.8, label.offset = 0.5, edge.width = 1.5)
tiplabels(pie = nest.mat, piecol = nest.cols, cex = 0.35)
nodelabels(pie = anc.states.m9, piecol = nest.cols, cex = 0.4)
legend(x = -2.25, y = 0, xjust = 0, yjust = 0,
       title = "NEST TYPE", title.adj = 0,
       legend = c("excavates burrow","adopts burrow or cavity",
                  "mud nest: open cup","mud nest: enclosed cup","mud nest: retort"),
       pch = 21, pt.bg = c("#2c7bb6","#abd9e9","#ffffbf","#fdae61","#d7191c"),
       pt.cex = 1.5, bty = "n", cex = 0.8)


library(rgdal)
library(rgeos)
library(rnaturalearth)
library(rworldmap)
library(maptools)
library(raster)
library(ape)
library(pryr)

wd = np("~/Google Drive/My Drive/projects/hirundinidae_phylogeny/workflow_reanalysis/analysis/biogeography/")
setwd(wd)
getwd()

#### Color palette ####
red <- "#d53e4f"
orange <- "#f46d43"
ltorange <- "#fdae61"
yellow <- "#fee08b"
grnyellow <- "#e6f598"
ltgreen <- "#abdda4"
green <- "#66c2a5"
blue <- "#3288bd"


#### Get some basemaps ####
land110 <- ne_download(scale = 110, type = 'land', category = 'physical')
lakes110 <- ne_download(scale = 110, type = 'lakes', category = 'physical')
oceans110 <- ne_download(scale = 110, type = 'ocean', category = 'physical')
countries110 <- ne_countries()
countries110 <- countries110[countries110$sovereignt != "Antarctica",]
land110 <- gIntersection(land110, countries110)

data("countriesLow"); map <- countriesLow
colSums(is.na(map@data))
map <- map[!is.na(map$continent),]
map <- map[!is.na(map$NAME),]
map@proj4string <- land110@proj4string


#### CMEC realms ####
cmec <- readOGR(dsn = "CMEC_regions_realms", layer = "newRealms")
levels(cmec@data$Realm)
# [1] "Afrotropical" "Australian" "Madagascan" "Nearctic" "Neotropical" "Oceanina"      
# [7] "Oriental" "Palearctic" "Panamanian" "Saharo-Arabian" "Sino-Japanese"

cmec.orig <- cmec

# Fix "Oceanina" --> "Oceanian
levels(cmec$Realm) <- c(levels(cmec$Realm), "Oceanian")
cmec$Realm[cmec$Realm == "Oceanina"] <- "Oceanian"
cmec@data <- droplevels(cmec@data)

# "Saharo-Arabian","Sino-Japanese", "Nearctic" --> "Palearctic"
# "Madagascan" --> "Afrotropical"
# "Panamanian" --> "Neotropical"
cmec$Realm[cmec$Realm == "Saharo-Arabian"] <- "Palearctic"
cmec$Realm[cmec$Realm == "Sino-Japanese"] <- "Palearctic"
cmec$Realm[cmec$Realm == "Nearctic"] <- "Palearctic"
cmec$Realm[cmec$Realm == "Madagascan"] <- "Afrotropical"
cmec$Realm[cmec$Realm == "Panamanian"] <- "Neotropical"
cmec@data <- droplevels(cmec@data)
levels(cmec$Realm)

# Merge the Palearctic, Neotropical, and Afrotropical polygons
temp <- gUnaryUnion(cmec, id = cmec@data$Realm)

# To add the data frame back, the row names must match
row.names(temp) = as.character(1:length(temp))

# Get the data to add back to the merged spatial polygons object
dat = unique(cmec@data$Realm)
dat = as.data.frame(dat)
colnames(dat) = "Realm"
dat$Realm <- sort(dat$Realm) # get it in a coherent order

# And add the data back to the map
temp = SpatialPolygonsDataFrame(temp, dat)

# Save the new object
cmec <- temp


## Pull out the realms as objects and clean them up
cmec@data

oldworld <- map[map$continent == "Africa" | map$continent == "Eurasia" | map$continent == "Australia",]
oldworld <- oldworld[oldworld$NAME != "Greenland",]

# Split the palearctic and nearctic, add greenland to the nearctic
temp <- cmec[cmec$Realm == "Palearctic",]
temp <- gIntersection(temp, land110)

nearctic <- gIntersection(temp, map[map$continent == "North America" | map$GEO3 == "Meso-America",])
nearctic <- spRbind(nearctic, countries110[countries110$name == "Greenland",])
plot(countries110, col = "lightblue", border = "lightblue")
plot(nearctic, add = T, lwd = 0.1, border = "lightgray", col = "#f46d43")

palearctic <- gIntersection(temp, oldworld)
plot(countries110, col = "lightblue", border = "lightblue")
plot(palearctic, add = T, lwd = 0.1, border = "lightgray", col = "#f46d43")

# Split the neotropical and panamanian
temp <- cmec[cmec$Realm == "Neotropical",]
temp <- gIntersection(temp, land110)

neotropical <- gIntersection(temp, map[map$GEO3 == "South America",])
plot(countries110, col = "lightblue", border = "lightblue")
plot(neotropical, add = T, lwd = 0.1, border = "lightgray", col = "#f46d43")

panamanian <- gIntersection(temp, map[map$GEO3 == "Caribbean" | map$GEO3 == "Meso-America",])
plot(countries110, col = "lightblue", border = "lightblue")
plot(panamanian, add = T, lwd = 0.1, border = "lightgray", col = "#f46d43")

# Afrotropical
afrotropical <- cmec[cmec$Realm == "Afrotropical",]
afrotropical <- gIntersection(land110, afrotropical)
afrotropical <- gIntersection(afrotropical, map[map$continent == "Africa",])
plot(countries110, col = "lightblue", border = "lightblue")
plot(afrotropical, add = T, lwd = 0.1, border = "lightgray", col = "#f46d43")

# Australian
australian <- cmec[cmec$Realm == "Australian",]
australian <- gIntersection(land110, australian)
plot(countries110, col = "lightblue", border = "lightblue")
plot(australian, add = T, lwd = 0.1, border = "lightgray", col = "#f46d43")

# Oceanian
oceanian <- cmec[cmec$Realm == "Oceanian",]
oceanian <- gIntersection(land110, oceanian)
plot(countries110, col = "lightblue", border = "lightblue")
plot(oceanian, add = T, lwd = 0.1, border = "lightgray", col = "#f46d43")

# Oriental
oriental <- cmec[cmec$Realm == "Oriental",]
oriental <- gIntersection(land110, oriental)
plot(countries110, col = "lightblue", border = "lightblue")
plot(oriental, add = T, lwd = 0.1, border = "lightgray", col = "#f46d43")

# Caspian Sea
caspian_bb <- as(extent(45, 56, 35, 50), "SpatialPolygons") # make a bounding box
proj4string(caspian_bb) <- countries110@proj4string
caspian <- gIntersection(oceans110, caspian_bb) # get the points that define the outline of the Caspian
plot(countries110, col = "lightblue", border = "lightblue")
plot(caspian, col = "red", border = "lightblue", lwd = 0.01, add = T)

# Plot everything together, make sure nothing weird is going on or missing.
# Color scheme
# Add letters

# axis(side = 1); axis(side = 2); grid()

realms.map %<a-% {
  plot(land110, border = "white", col = "lightgray")
  plot(caspian, col = "lightgray", border = "lightgray", lwd = 1, add = T)
  plot(nearctic, add = T, lwd = 0.25, border = green, col = green)
  plot(palearctic, add = T, lwd = 0.25, border = blue, col = blue)
  plot(neotropical, add = T, lwd = 0.25, border = orange, col = orange)
  plot(panamanian, add = T, lwd = 0.1, border = grnyellow, col = grnyellow)
  plot(afrotropical, add = T, lwd = 0.1, border = yellow, col = yellow)
  plot(australian, add = T, lwd = 0.1, border = ltorange, col = ltorange)
  plot(oceanian, add = T, lwd = 0.1, border = red, col = red)
  plot(oriental, add = T, lwd = 0.1, border = ltgreen, col = ltgreen)
  plot(land110, border = "gray33", add = T)
  text("A", x=20, y=10)
  text("D", x=135, y=-23)
  text("F", x=-100, y=45)
  text("H", x=-60, y=-10)
  text("E", x=145, y=5)
  text("C", x=78, y=24)
  text("B", x=85, y=55)
  text("G", x=-95, y=7)
}

realms.map


#### Plot the map of realms with the tree ####

tr <- read.tree("./hirundinidae-all-rev1-95-prune.tree")
as.matrix(tr$tip.label, ncol = 1)
#tr$tip.label <- c("Atticora_pileata","Atticora_tibialis","Atticora_fasciata","Pygochelidon_cyanoleuca","Pygochelidon_melanoleuca","Alopochelidon_fucata","Orochelidon_murina","Orochelidon_andecola","Orochelidon_flavipes","Stelgidopteryx_serripennis_ridgwayi","Stelgidopteryx_serripennis","Stelgidopteryx_ruficollis","Progne_tapera","Progne_subis","Progne_elegans","Progne_dominicensis","Progne_sinaloae","Progne_chalybea","Progne_murphyi","Tachycineta_bicolor","Tachycineta_cyaneoviridis","Tachycineta_euchrysea","Tachycineta_thalassina_brachyptera","Tachycineta_thalassina_thalassina","Tachycineta_leucorrhoa","Tachycineta_leucopyga","Tachycineta_stolzmanni","Tachycineta_albiventer","Tachycineta_albilinea","Riparia_riparia","Riparia_diluta","Riparia_paludicola_cowani","Riparia_paludicola_paludicola","Phedina_borbonica","Phedinopsis_brazzae","Neophedina_cincta","Cheramoeca_leucosterna","Pseudhirundo_griseopyga","Petrochelidon_fluvicola","Petrochelidon_ariel","Petrochelidon_nigricans","Petrochelidon_pyrrhonota","Petrochelidon_fulva","Petrochelidon_rufocollaris","Petrochelidon_preussi","Petrochelidon_spilodera","Cecropis_abyssinica","Cecropis_cucullata","Cecropis_daurica_Eurasia/SE_Asia","Cecropis_daurica_Africa","Cecropis_senegalensis","Cecropis_semirufa_semirufa","Cecropis_semirufa_gordoni","Petrochelidon_fuliginosa","Delichon_nipalense","Delichon_dasypus","Delichon_urbicum_lagopodum","Delichon_urbicum_urbicum","Hirundo_atrocaerulea","Hirundo_dimidiata","Hirundo_rustica","Hirundo_angolensis","Hirundo_lucida","Hirundo_aethiopica","Hirundo_albigularis","Hirundo_nigrita","Hirundo_smithii","Hirundo_neoxena","Hirundo_tahitica","Ptyonoprogne_fuligula_S._pop.","Ptyonoprogne_fuligula_N_pop.","Ptyonoprogne_concolor","Ptyonoprogne_rupestris","Psalidoprocne_pristoptera_N.","Psalidoprocne_albiceps","Psalidoprocne_pristoptera_S.","Psalidoprocne_fuliginosa","Psalidoprocne_obscura","Psalidoprocne_nitens","Pseudochelidon_sirintarae","Pseudochelidon_eurystomina")

pdf("./temp.pdf", height = 11, width = 8.5)

par(fig = c(0.05,1,0.2,1))
plot(tr, no.margin = T, cex = 0.6, label.offset = 0.7, edge.width = 1)
par(fig = c(0,1,0.025,0.225), new = T)
realms.map
dev.off()


#### Get the node and tip colors from the BioGeoBEARS DEC+j analysis ####

## This is a surprisingly irritating ordeal.
# Go through BioGeoBEARS_constrained_dispersal.R to the point where the saved resDECj object is loaded.
  # Line 464 in script.
summary(resDECj)

# And then through BioGeoBEARS_get_MLstates.R to get the ancestral states and the tip states
# This runs part of plot_BioGeoBEARS_results
# Output is a matrix of node and tip states, and a matrix of colors for those states.
MLstates
cols_byNode
Ntip(tr); Nnode(tr)

# Plot the states and colors on the tree
plot(tr, no.margin = T, cex = 0.6, label.offset = 0.5)
nodelabels(text = MLstates[89:169], bg = cols_byNode[89:169], cex = 0.6, font = 2)
tiplabels(text = MLstates[1:88], bg = cols_byNode[1:88], cex = 0.6, font = 2)


#### Plot the complete figure ####

#tr$tip.label <- c("Atticora_pileata","Atticora_tibialis","Atticora_fasciata","Pygochelidon_cyanoleuca","Pygochelidon_melanoleuca","Alopochelidon_fucata","Orochelidon_murina","Orochelidon_andecola","Orochelidon_flavipes","Stelgidopteryx_serripennis_ridgwayi","Stelgidopteryx_serripennis","Stelgidopteryx_ruficollis","Progne_tapera","Progne_subis","Progne_elegans","Progne_dominicensis","Progne_sinaloae","Progne_chalybea","Progne_murphyi","Tachycineta_bicolor","Tachycineta_cyaneoviridis","Tachycineta_euchrysea","Tachycineta_thalassina_brachyptera","Tachycineta_thalassina_thalassina","Tachycineta_leucorrhoa","Tachycineta_leucopyga","Tachycineta_stolzmanni","Tachycineta_albiventer","Tachycineta_albilinea","Riparia_riparia","Riparia_diluta","Riparia_paludicola_cowani","Riparia_paludicola_paludicola","Phedina_borbonica","Phedinopsis_brazzae","Neophedina_cincta","Cheramoeca_leucosterna","Pseudhirundo_griseopyga","Petrochelidon_fluvicola","Petrochelidon_ariel","Petrochelidon_nigricans","Petrochelidon_pyrrhonota","Petrochelidon_fulva","Petrochelidon_rufocollaris","Petrochelidon_preussi","Petrochelidon_spilodera","Cecropis_abyssinica","Cecropis_cucullata","Cecropis_daurica_Eurasia/SE_Asia","Cecropis_daurica_Africa","Cecropis_senegalensis","Cecropis_semirufa_semirufa","Cecropis_semirufa_gordoni","Petrochelidon_fuliginosa","Delichon_nipalense","Delichon_dasypus","Delichon_urbicum_lagopodum","Delichon_urbicum_urbicum","Hirundo_atrocaerulea","Hirundo_dimidiata","Hirundo_rustica","Hirundo_angolensis","Hirundo_lucida","Hirundo_aethiopica","Hirundo_albigularis","Hirundo_nigrita","Hirundo_smithii","Hirundo_neoxena","Hirundo_tahitica","Ptyonoprogne_fuligula_S._pop.","Ptyonoprogne_fuligula_N_pop.","Ptyonoprogne_concolor","Ptyonoprogne_rupestris","Psalidoprocne_pristoptera_N.","Psalidoprocne_albiceps","Psalidoprocne_pristoptera_S.","Psalidoprocne_fuliginosa","Psalidoprocne_obscura","Psalidoprocne_nitens","Pseudochelidon_sirintarae","Pseudochelidon_eurystomina")


# Map at the bottom, no time axis
pdf("./temp.pdf", height = 11, width = 8.5)
par(fig = c(0.05,1,0.2,1))
plot(tr, no.margin = T, cex = 0.6, label.offset = 1.0)
nodelabels(text = MLstates[89:169], bg = cols_byNode[89:169], cex = 0.6, font = 2)
tiplabels(text = MLstates[1:88], bg = cols_byNode[1:88], cex = 0.6, font = 2)
par(fig = c(0,1,0.025,0.225), new = T)
realms.map
dev.off()

# Map at the top, with time axis
pdf("./temp.pdf", height = 11, width = 8.5)
par(fig = c(0.05,1,0,0.8), mar = c(2,0,0,0))
plot(tr, no.margin = F, cex = 0.6, label.offset = 0.6, edge.width = 1.2)
nodelabels(text = MLstates[89:169], bg = cols_byNode[89:169], cex = 0.5, font = 2)
tiplabels(text = MLstates[1:88], bg = cols_byNode[1:88], cex = 0.5, font = 2)
axis(side = 1, lwd = 1.3, tcl = 0.5, cex.axis = 0.8, font = 2, padj = -2, pos = -1.5,
     at = c(0, 5, 10, 15), labels = c(15, 10, 5, 0))
mtext("Ma", at = 15.3, side = 1, cex = 0.8, padj = -2.25, font = 2, adj = 0)
par(fig = c(0,0.75,0.745,0.975), new = T)
realms.map
dev.off()

# Map at the top, with time axis, arrows for cavity adoption and mud-nesting nodes
pdf("../../plots/biogeography_DIVALIKE+j_map_scale.pdf", height = 14, width = 10)
par(fig = c(0.05,1,0,0.8), mar = c(2,0,0,0))
plot(tr, no.margin = F, cex = 0.6, label.offset = 0.6, edge.width = 1.2)
nodelabels(text = MLstates[89:169], bg = cols_byNode[89:169], cex = 0.5, font = 2)
tiplabels(text = MLstates[1:88], bg = cols_byNode[1:88], cex = 0.5, font = 2)
axis(side = 1, lwd = 1.3, tcl = 0.5, cex.axis = 0.8, font = 2, padj = -2, pos = -1.5,
     at = c(0, 5, 10, 15), labels = c(15, 10, 5, 0))
mtext("Ma", at = 15.3, side = 1, cex = 0.8, padj = -2.25, font = 2, adj = 0)
arrows(3.25,62,4.3,59.5, col = "red", lwd = 2, angle = 20, length = 0.15)
arrows(4.8,13,5.85,15.75, col = "blue", lwd = 2, angle = 20, length = 0.15)

par(fig = c(0,0.75,0.745,0.975), new = T)
realms.map
dev.off()

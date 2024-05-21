library(ape)
library(phytools)

outgroup <- c("Erythrocercus_mccallii_KU_8688", "Locustella_lanceolata_KU_4248")

### Concatenated maximum likelihood tree (ha95-rev1 dataset; n = 118 taxa)

# Read in and root tree; prune outgroups
ha95.rax <- read.tree("~/Google Drive/My Drive/projects/hirundinidae_phylogeny/workflow_reanalysis/analysis/concatenated/ha95-rev1/ha95-rev1-annotated.raxml.support.tree")
#write(ha95.rax$tip.label,'~/Google Drive/My Drive/projects/hirundinidae_phylogeny/workflow_reanalysis/plots/ha95-raxml-tips.txt') # Edit separate taxon and locality/clade labels from these below
ha95.rax <- root(ha95.rax, outgroup, edgelabel = T)
keep <- c("Pseudochelidon_eurystomina_Congo_FMNH_213526",
          "Psalidoprocne_albiceps_Malawi_FMNH_467949",
          "Ptyonoprogne_rupestris_Greece_LSUMNS_25377",
          "Hirundo_smithii_Ghana_LSUMNS_39509",
          "Petrochelidon_fuliginosa_Cameroon_KU_131409",
          "Delichon_urbicum_urbicum_South_Africa_LSUMNS_14053",
          "Petrochelidon_pyrrhonota_USA_LSUMNS_63941",
          "Cecropis_cucullata_South_Africa_LSUMNS_14063",
          "Cheramoeca_leucosterna_Australia_LSUMNS_14165",
          "Pseudhirundo_griseopyga_South_Africa_LSUMNS_34230",
          "Phedina_borbonica_Madagascar_LSUMNS_25388",
          "Riparia_riparia_eastern_Palearctic_Russia_UWBM_46946",
          "Tachycineta_thalassina_thalassina_USA_MVZ_182222",
          "Stelgidopteryx_ruficollis_Peru_Peru_LSUMNS_74771",
          "Atticora_fasciata_Peru_LSUMNS_75904",
          "Pygochelidon_cyanoleuca_Peru_LSUMNS_43665",
          "Alopochelidon_fucata_Uruguay_CUMV_50652",
          "Orochelidon_murina_Peru_LSUMNS_32352",
          "Progne_tapera_tapera_Peru_LSUMNS_75925")

ha95.rax <- keep.tip(ha95.rax, keep)
ha95.rax$edge.length <- NULL
ha95.rax <- ladderize(ha95.rax)
plot(ha95.rax, cex = 0.7, no.margin = T, label.offset = 0.15,edge.width = 1)

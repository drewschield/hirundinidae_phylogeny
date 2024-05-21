library(treeio)
#library(tracerer)
#library(ape)
#library(phytools)
#library(devtools)

# We're going to use treeio (https://yulab-smu.top/treedata-book/chapter2.html) to import, prune, and export the time-calibrated phylogeny.

### Bayesian time-calibrated phylogeny (ha95-rev1 dataset; n = 118 taxa)

# Read in and root tree; prune to 88 representative taxa
ha95.beast <- read.beast("~/Google Drive/My Drive/projects/hirundinidae_phylogeny/workflow_reanalysis/analysis/beast/hirundinidae-all-rev1-95.tree")
prune <- c("Erythrocercus_mccallii_KU_8688",
           "Locustella_lanceolata_KU_4248",
           "Hirundo_rustica_erythrogaster_USA_LSUMNS_20587",
           "Hirundo_rustica_gutturalis_Russia_UWBM_47365",
           "Hirundo_rustica_tytleri_Russia_UWBM_51725",
           "Petrochelidon_fulva_Caribbean_clade_Jamaica_LSUMNS_27022",
           "Cecropis_abyssinica_unitatis_South_Africa_LSUMNS_34246",
           "Cecropis_daurica_Kuwait_Kuwait_LSUMNS_87323",
           "Cecropis_daurica_Russia_Russia_UWBM_75471",
           "Cecropis_daurica_Greece_Greece_YPM_ORN_142643",
           "Cecropis_daurica_Pakistan_Pakistan_FMNH_347828",
           "Cecropis_daurica_Singapore_Singapore_UWBM_83563",
           "Cecropis_striolata_Philippines_Philippines_KU_119813",
           "Cecropis_striolata_Vietnam_Vietnam_KU_119705",
           "Riparia_diluta_Russia_Russia_UWBM_67556",
           "Riparia_riparia_North_America_USA_LSUMNS_33324",
           "Riparia_riparia_western_Palearctic_Russia_UWBM_82230",
           "Stelgidopteryx_ruficollis_decolor_Costa_Rica_LSUMNS_27259",
           "Stelgidopteryx_serripennis_northern_clade_USA_LSUMNS_57224",
           "Stelgidopteryx_serripennis_ridgwayi_Mexico_LSUMNS_27213",
           "Progne_chalybea_Amazon_Bolivia_LSUMNS_79812",
           "Progne_chalybea_Amazon_Brazil_LSUMNS_86649",
           "Progne_chalybea_Central_America_Panama_LSUMNS_28811",
           "Progne_chalybea_Tumbes_Peru_LSUMNS_66042",
           "Progne_dominicensis_Dominican_Republic_LSUMNS_22020",
           "Progne_elegans_Brazil_LSUMNS_25506",
           "Progne_murphyi_Peru_LSUMZ_114186",
           "Progne_sinaloae_Mexico_KU_40045",
           "Progne_subis_USA_LSUMNS_41548",
           "Progne_tapera_fusca_Bolivia_LSUMNS_37913")

ha95.beast <- drop.tip(ha95.beast,tip = prune)
write.beast(ha95.beast,file='~/Google Drive/My Drive/projects/hirundinidae_phylogeny/workflow_reanalysis/analysis/beast/hirundinidae-all-rev1-95-prune.tree')


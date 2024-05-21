# Hirundinidae Phylogeny

![Consensus tree](/tree.png "cover image")

This repository contains details on the data processing and analysis steps used to estimate the swallow family phylogeny, divergence times, historical biogeography, and ancestral state reconstructions for sociality and nesting type. This workflow is a companion to the methods described in Schield & Brown et al. (_in Press_).

Much of the code and scripts here were adopted from elements of co-author Clare Brown's notes (thank you, Clare!).

Note: this repository assumes a specific file organization. Your environment may vary and will need to be adjusted accordingly.

## Part 0 - General information and installation of phyluce

The dataset for this study consists of UCE loci sampled from a majority of swallow species.

The samples include preserved tissue and historical skins, the latter which yielded lower quality genomic data, resulting in spurious results and very long branch lengths for some samples from skins (and several that were dropped from analysis completely).

------------------------------------------------------------------------------------------
### Set up environment

```
cd /data3/hirundinidae_phylogeny
mkdir workflow_reanalysis
```

------------------------------------------------------------------------------------------
### Install phyluce

phyluce is the software package and wrapper that we'll use to process, assemble, correct, and align UCE loci.

It is installed via Conda; Install it in it's own environment from the Linux release version on the phyluce GitHub page (https://github.com/faircloth-lab/phyluce/releases).
```
cd ~/tmp/
mkdir phyluce
cd phyluce
wget https://raw.githubusercontent.com/faircloth-lab/phyluce/v1.7.3/distrib/phyluce-1.7.3-py36-Linux-conda.yml
conda env create -n phyluce --file phyluce-1.7.3-py36-Linux-conda.yml
```

To activate/deactivate environment:
```
conda activate phyluce
conda deactivate
```

------------------------------------------------------------------------------------------
### Install CIAlign

Use CIAlign to diagnose issues with skin samples in downstream alignments.

```
conda create -n "cialign" python=3.6
conda activate cialign
conda install -c bioconda cialign
```

To activate/deactivate environment:
```
conda activate cialign
conda deactivate
```

## Part 1 - Processing and correction of UCE data for skin samples

Revisit the input data for skin samples to determine if running the new correction step in phyluce improves things, and to see if we can include a greater number of these samples in the final analyses.

Analysis steps that have already been completed:
* Initial processing (cleaning and quantification of raw read data)
* SPAdes assembly of UCE contigs
* Find UCEs

Note: sample Progne_modesta_Ecuador_MVZ_130128 has already been completely removed from analysis at this point due to a very small number of assembled contigs.

Here, we will pick up with the assembled contigs for the remaining 22 skin samples, and run mapping and correction in phyluce using the contigs and fastq files for each.

The skin samples are:

* Cecropis_daurica_Africa_Uganda_LACM_71522
* Delichon_nipalense_Nepal_FMNH_276821
* Hirundo_lucida_Gambia_FMNH_320046
* Hirundo_megaensis_Ethiopia_AMNH_SKIN_348889
* Hirundo_nigrorufa_Angola_AMNH_SKIN_707943
* Petrochelidon_fluvicola_India_FMNH_233414
* Petrochelidon_rufigula_Angola_AMNH_SKIN_707945
* Phedinopsis_brazzae_Angola_AMNH_SKIN_764767
* Progne_cryptoleuca_Cuba_LSUMZ_142940
* Progne_cryptoleuca_Cuba_LSUMZ_142944
* Progne_modesta_Ecuador_MVZ_130129
* Progne_murphyi_Peru_LSUMZ_114185
* Progne_murphyi_Peru_LSUMZ_114186
* Progne_sinaloae_Mexico_KU_40044
* Progne_sinaloae_Mexico_KU_40045
* Psalidoprocne_obscura_Ivory_Coast_FMNH_279082
* Pseudochelidon_eurystomina_Congo_FMNH_213526
* Pseudochelidon_sirintarae_Thailand_AMNH_SKIN_708673
* Ptyonoprogne_concolor_India_FMNH_233366
* Ptyonoprogne_obsoleta_Egypt_UMMZ_224076
* Riparia_congica_Congo_FMNH_213543

------------------------------------------------------------------------------------------
### Input data locations

Contig fasta files: `/data3/hirundinidae_phylogeny/Hirundinidae_UCE/phyluce/spades-assemblies-ingroup-skins/contigs`. Note: the .fasta files here are aliases to contigs.fasta files in respective sample-specific directories (in the same parent directory as `.../contigs/`)

Fastq files: `/data3/hirundinidae_phylogeny/Hirundinidae_UCE/phyluce/clean-fastq-ingroup/<sample>/split-adapter-quality-trimmed`

------------------------------------------------------------------------------------------
### Set up environment

We'll make a separate `phyluce` directory in `/data3/hirundinidae_phylogeny/workflow_reanalysis/` for mapping and correction of the 21 skin samples.

```
cd /data3/hirundinidae_phylogeny/workflow_reanalysis/
mkdir phyluce
cd phyluce
```

------------------------------------------------------------------------------------------
### 1. Mapping fastq data to contigs using phyluce

We'll use the [mapping](https://phyluce.readthedocs.io/en/latest/daily-use/daily-use-4-workflows.html#mapping) procedure in phyluce to generate .bam files for each of the skin samples.

This requires an input config YAML file like this example:
```
reads:
    alligator-mississippiensis: ../../phyluce/tests/test-expected/raw-reads/alligator-mississippiensis/
    gallus-gallus: ../../phyluce/tests/test-expected/raw-reads/gallus-gallus
    peromyscus-maniculatus: ../../phyluce/tests/test-expected/raw-reads/peromyscus-maniculatus
    rana-sphenocephafa: ../../phyluce/tests/test-expected/raw-reads/rana-sphenocephafa

contigs:
    alligator-mississippiensis: ../../phyluce/tests/test-expected/spades/contigs/alligator_mississippiensis.contigs.fasta
    gallus-gallus: ../../phyluce/tests/test-expected/spades/contigs/gallus_gallus.contigs.fasta
    peromyscus-maniculatus: ../../phyluce/tests/test-expected/spades/contigs/peromyscus_maniculatus.contigs.fasta
    rana-sphenocephafa: ../../phyluce/tests/test-expected/spades/contigs/rana_sphenocephafa.contigs.fasta
```

#### Format config YAML for skin samples

Note: do not use TAB indentation!

1-mapping.config
```
reads:
    Cecropis_daurica_Africa_Uganda_LACM_71522: /data3/hirundinidae_phylogeny/Hirundinidae_UCE/phyluce/clean-fastq-ingroup/Cecropis_daurica_Africa_Uganda_LACM_71522/split-adapter-quality-trimmed
    Delichon_nipalense_Nepal_FMNH_276821: /data3/hirundinidae_phylogeny/Hirundinidae_UCE/phyluce/clean-fastq-ingroup/Delichon_nipalense_Nepal_FMNH_276821/split-adapter-quality-trimmed
    Hirundo_lucida_Gambia_FMNH_320046: /data3/hirundinidae_phylogeny/Hirundinidae_UCE/phyluce/clean-fastq-ingroup/Hirundo_lucida_Gambia_FMNH_320046/split-adapter-quality-trimmed
    Hirundo_megaensis_Ethiopia_AMNH_SKIN_348889: /data3/hirundinidae_phylogeny/Hirundinidae_UCE/phyluce/clean-fastq-ingroup/Hirundo_megaensis_Ethiopia_AMNH_SKIN_348889/split-adapter-quality-trimmed
    Hirundo_nigrorufa_Angola_AMNH_SKIN_707943: /data3/hirundinidae_phylogeny/Hirundinidae_UCE/phyluce/clean-fastq-ingroup/Hirundo_nigrorufa_Angola_AMNH_SKIN_707943/split-adapter-quality-trimmed
    Petrochelidon_fluvicola_India_FMNH_233414: /data3/hirundinidae_phylogeny/Hirundinidae_UCE/phyluce/clean-fastq-ingroup/Petrochelidon_fluvicola_India_FMNH_233414/split-adapter-quality-trimmed
    Petrochelidon_rufigula_Angola_AMNH_SKIN_707945: /data3/hirundinidae_phylogeny/Hirundinidae_UCE/phyluce/clean-fastq-ingroup/Petrochelidon_rufigula_Angola_AMNH_SKIN_707945/split-adapter-quality-trimmed
    Phedinopsis_brazzae_Angola_AMNH_SKIN_764767: /data3/hirundinidae_phylogeny/Hirundinidae_UCE/phyluce/clean-fastq-ingroup/Phedinopsis_brazzae_Angola_AMNH_SKIN_764767/split-adapter-quality-trimmed
    Progne_cryptoleuca_Cuba_LSUMZ_142940: /data3/hirundinidae_phylogeny/Hirundinidae_UCE/phyluce/clean-fastq-ingroup/Progne_cryptoleuca_Cuba_LSUMZ_142940/split-adapter-quality-trimmed
    Progne_cryptoleuca_Cuba_LSUMZ_142944: /data3/hirundinidae_phylogeny/Hirundinidae_UCE/phyluce/clean-fastq-ingroup/Progne_cryptoleuca_Cuba_LSUMZ_142944/split-adapter-quality-trimmed
    Progne_modesta_Ecuador_MVZ_130129: /data3/hirundinidae_phylogeny/Hirundinidae_UCE/phyluce/clean-fastq-ingroup/Progne_modesta_Ecuador_MVZ_130129/split-adapter-quality-trimmed
    Progne_murphyi_Peru_LSUMZ_114185: /data3/hirundinidae_phylogeny/Hirundinidae_UCE/phyluce/clean-fastq-ingroup/Progne_murphyi_Peru_LSUMZ_114185/split-adapter-quality-trimmed
    Progne_murphyi_Peru_LSUMZ_114186: /data3/hirundinidae_phylogeny/Hirundinidae_UCE/phyluce/clean-fastq-ingroup/Progne_murphyi_Peru_LSUMZ_114186/split-adapter-quality-trimmed
    Progne_sinaloae_Mexico_KU_40044: /data3/hirundinidae_phylogeny/Hirundinidae_UCE/phyluce/clean-fastq-ingroup/Progne_sinaloae_Mexico_KU_40044/split-adapter-quality-trimmed
    Progne_sinaloae_Mexico_KU_40045: /data3/hirundinidae_phylogeny/Hirundinidae_UCE/phyluce/clean-fastq-ingroup/Progne_sinaloae_Mexico_KU_40045/split-adapter-quality-trimmed
    Psalidoprocne_obscura_Ivory_Coast_FMNH_279082: /data3/hirundinidae_phylogeny/Hirundinidae_UCE/phyluce/clean-fastq-ingroup/Psalidoprocne_obscura_Ivory_Coast_FMNH_279082/split-adapter-quality-trimmed
    Pseudochelidon_eurystomina_Congo_FMNH_213526: /data3/hirundinidae_phylogeny/Hirundinidae_UCE/phyluce/clean-fastq-ingroup/Pseudochelidon_eurystomina_Congo_FMNH_213526/split-adapter-quality-trimmed
    Pseudochelidon_sirintarae_Thailand_AMNH_SKIN_708673: /data3/hirundinidae_phylogeny/Hirundinidae_UCE/phyluce/clean-fastq-ingroup/Pseudochelidon_sirintarae_Thailand_AMNH_SKIN_708673/split-adapter-quality-trimmed
    Ptyonoprogne_concolor_India_FMNH_233366: /data3/hirundinidae_phylogeny/Hirundinidae_UCE/phyluce/clean-fastq-ingroup/Ptyonoprogne_concolor_India_FMNH_233366/split-adapter-quality-trimmed
    Ptyonoprogne_obsoleta_Egypt_UMMZ_224076: /data3/hirundinidae_phylogeny/Hirundinidae_UCE/phyluce/clean-fastq-ingroup/Ptyonoprogne_obsoleta_Egypt_UMMZ_224076/split-adapter-quality-trimmed
    Riparia_congica_Congo_FMNH_213543: /data3/hirundinidae_phylogeny/Hirundinidae_UCE/phyluce/clean-fastq-ingroup/Riparia_congica_Congo_FMNH_213543/split-adapter-quality-trimmed
contigs:
    Cecropis_daurica_Africa_Uganda_LACM_71522: /data3/hirundinidae_phylogeny/Hirundinidae_UCE/phyluce/spades-assemblies-ingroup-skins/Cecropis_daurica_Africa_Uganda_LACM_71522_spades/contigs.fasta
    Delichon_nipalense_Nepal_FMNH_276821: /data3/hirundinidae_phylogeny/Hirundinidae_UCE/phyluce/spades-assemblies-ingroup-skins/Delichon_nipalense_Nepal_FMNH_276821_spades/contigs.fasta
    Hirundo_lucida_Gambia_FMNH_320046: /data3/hirundinidae_phylogeny/Hirundinidae_UCE/phyluce/spades-assemblies-ingroup-skins/Hirundo_lucida_Gambia_FMNH_320046_spades/contigs.fasta
    Hirundo_megaensis_Ethiopia_AMNH_SKIN_348889: /data3/hirundinidae_phylogeny/Hirundinidae_UCE/phyluce/spades-assemblies-ingroup-skins/Hirundo_megaensis_Ethiopia_AMNH_SKIN_348889_spades/contigs.fasta
    Hirundo_nigrorufa_Angola_AMNH_SKIN_707943: /data3/hirundinidae_phylogeny/Hirundinidae_UCE/phyluce/spades-assemblies-ingroup-skins/Hirundo_nigrorufa_Angola_AMNH_SKIN_707943_spades/contigs.fasta
    Petrochelidon_fluvicola_India_FMNH_233414: /data3/hirundinidae_phylogeny/Hirundinidae_UCE/phyluce/spades-assemblies-ingroup-skins/Petrochelidon_fluvicola_India_FMNH_233414_spades/contigs.fasta
    Petrochelidon_rufigula_Angola_AMNH_SKIN_707945: /data3/hirundinidae_phylogeny/Hirundinidae_UCE/phyluce/spades-assemblies-ingroup-skins/Petrochelidon_rufigula_Angola_AMNH_SKIN_707945_spades/contigs.fasta
    Phedinopsis_brazzae_Angola_AMNH_SKIN_764767: /data3/hirundinidae_phylogeny/Hirundinidae_UCE/phyluce/spades-assemblies-ingroup-skins/Phedinopsis_brazzae_Angola_AMNH_SKIN_764767_spades/contigs.fasta
    Progne_cryptoleuca_Cuba_LSUMZ_142940: /data3/hirundinidae_phylogeny/Hirundinidae_UCE/phyluce/spades-assemblies-ingroup-skins/Progne_cryptoleuca_Cuba_LSUMZ_142940_spades/contigs.fasta
    Progne_cryptoleuca_Cuba_LSUMZ_142944: /data3/hirundinidae_phylogeny/Hirundinidae_UCE/phyluce/spades-assemblies-ingroup-skins/Progne_cryptoleuca_Cuba_LSUMZ_142944_spades/contigs.fasta
    Progne_modesta_Ecuador_MVZ_130129: /data3/hirundinidae_phylogeny/Hirundinidae_UCE/phyluce/spades-assemblies-ingroup-skins/Progne_modesta_Ecuador_MVZ_130129_spades/contigs.fasta
    Progne_murphyi_Peru_LSUMZ_114185: /data3/hirundinidae_phylogeny/Hirundinidae_UCE/phyluce/spades-assemblies-ingroup-skins/Progne_murphyi_Peru_LSUMZ_114185_spades/contigs.fasta
    Progne_murphyi_Peru_LSUMZ_114186: /data3/hirundinidae_phylogeny/Hirundinidae_UCE/phyluce/spades-assemblies-ingroup-skins/Progne_murphyi_Peru_LSUMZ_114186_spades/contigs.fasta
    Progne_sinaloae_Mexico_KU_40044: /data3/hirundinidae_phylogeny/Hirundinidae_UCE/phyluce/spades-assemblies-ingroup-skins/Progne_sinaloae_Mexico_KU_40044_spades/contigs.fasta
    Progne_sinaloae_Mexico_KU_40045: /data3/hirundinidae_phylogeny/Hirundinidae_UCE/phyluce/spades-assemblies-ingroup-skins/Progne_sinaloae_Mexico_KU_40045_spades/contigs.fasta
    Psalidoprocne_obscura_Ivory_Coast_FMNH_279082: /data3/hirundinidae_phylogeny/Hirundinidae_UCE/phyluce/spades-assemblies-ingroup-skins/Psalidoprocne_obscura_Ivory_Coast_FMNH_279082_spades/contigs.fasta
    Pseudochelidon_eurystomina_Congo_FMNH_213526: /data3/hirundinidae_phylogeny/Hirundinidae_UCE/phyluce/spades-assemblies-ingroup-skins/Pseudochelidon_eurystomina_Congo_FMNH_213526_spades/contigs.fasta
    Pseudochelidon_sirintarae_Thailand_AMNH_SKIN_708673: /data3/hirundinidae_phylogeny/Hirundinidae_UCE/phyluce/spades-assemblies-ingroup-skins/Pseudochelidon_sirintarae_Thailand_AMNH_SKIN_708673_spades/contigs.fasta
    Ptyonoprogne_concolor_India_FMNH_233366: /data3/hirundinidae_phylogeny/Hirundinidae_UCE/phyluce/spades-assemblies-ingroup-skins/Ptyonoprogne_concolor_India_FMNH_233366_spades/contigs.fasta
    Ptyonoprogne_obsoleta_Egypt_UMMZ_224076: /data3/hirundinidae_phylogeny/Hirundinidae_UCE/phyluce/spades-assemblies-ingroup-skins/Ptyonoprogne_obsoleta_Egypt_UMMZ_224076_spades/contigs.fasta
    Riparia_congica_Congo_FMNH_213543: /data3/hirundinidae_phylogeny/Hirundinidae_UCE/phyluce/spades-assemblies-ingroup-skins/Riparia_congica_Congo_FMNH_213543_spades/contigs.fasta
```

#### Run phyluce mapping workflow

```
conda activate phyluce
phyluce_workflow --config ./1-mapping.config --output ./1-mapping --workflow mapping --cores 16
```

#### Output files

The outputs from the mapping step are in:
```
./1-mapping/coverage
./1-mapping/mapped_reads
./1-mapping/references
```

------------------------------------------------------------------------------------------
### 2. Correction of low depth, low quality bases in contigs

We'll use the [correction](https://phyluce.readthedocs.io/en/latest/daily-use/daily-use-4-workflows.html#correction) procedure in phyluce to generate fixed contig fasta files for each of the skin samples.

#### Format YAML config file with .bam and contig file paths

2-correction.config
```
bams:
    Cecropis_daurica_Africa_Uganda_LACM_71522: /data3/hirundinidae_phylogeny/workflow_reanalysis/phyluce/1-mapping/mapped_reads/Cecropis_daurica_Africa_Uganda_LACM_71522.fxm.sorted.md.bam 
    Delichon_nipalense_Nepal_FMNH_276821: /data3/hirundinidae_phylogeny/workflow_reanalysis/phyluce/1-mapping/mapped_reads/Delichon_nipalense_Nepal_FMNH_276821.fxm.sorted.md.bam 
    Hirundo_lucida_Gambia_FMNH_320046: /data3/hirundinidae_phylogeny/workflow_reanalysis/phyluce/1-mapping/mapped_reads/Hirundo_lucida_Gambia_FMNH_320046.fxm.sorted.md.bam 
    Hirundo_megaensis_Ethiopia_AMNH_SKIN_348889: /data3/hirundinidae_phylogeny/workflow_reanalysis/phyluce/1-mapping/mapped_reads/Hirundo_megaensis_Ethiopia_AMNH_SKIN_348889.fxm.sorted.md.bam 
    Hirundo_nigrorufa_Angola_AMNH_SKIN_707943: /data3/hirundinidae_phylogeny/workflow_reanalysis/phyluce/1-mapping/mapped_reads/Hirundo_nigrorufa_Angola_AMNH_SKIN_707943.fxm.sorted.md.bam 
    Petrochelidon_fluvicola_India_FMNH_233414: /data3/hirundinidae_phylogeny/workflow_reanalysis/phyluce/1-mapping/mapped_reads/Petrochelidon_fluvicola_India_FMNH_233414.fxm.sorted.md.bam 
    Petrochelidon_rufigula_Angola_AMNH_SKIN_707945: /data3/hirundinidae_phylogeny/workflow_reanalysis/phyluce/1-mapping/mapped_reads/Petrochelidon_rufigula_Angola_AMNH_SKIN_707945.fxm.sorted.md.bam 
    Phedinopsis_brazzae_Angola_AMNH_SKIN_764767: /data3/hirundinidae_phylogeny/workflow_reanalysis/phyluce/1-mapping/mapped_reads/Phedinopsis_brazzae_Angola_AMNH_SKIN_764767.fxm.sorted.md.bam 
    Progne_cryptoleuca_Cuba_LSUMZ_142940: /data3/hirundinidae_phylogeny/workflow_reanalysis/phyluce/1-mapping/mapped_reads/Progne_cryptoleuca_Cuba_LSUMZ_142940.fxm.sorted.md.bam 
    Progne_cryptoleuca_Cuba_LSUMZ_142944: /data3/hirundinidae_phylogeny/workflow_reanalysis/phyluce/1-mapping/mapped_reads/Progne_cryptoleuca_Cuba_LSUMZ_142944.fxm.sorted.md.bam 
    Progne_modesta_Ecuador_MVZ_130129: /data3/hirundinidae_phylogeny/workflow_reanalysis/phyluce/1-mapping/mapped_reads/Progne_modesta_Ecuador_MVZ_130129.fxm.sorted.md.bam 
    Progne_murphyi_Peru_LSUMZ_114185: /data3/hirundinidae_phylogeny/workflow_reanalysis/phyluce/1-mapping/mapped_reads/Progne_murphyi_Peru_LSUMZ_114185.fxm.sorted.md.bam 
    Progne_murphyi_Peru_LSUMZ_114186: /data3/hirundinidae_phylogeny/workflow_reanalysis/phyluce/1-mapping/mapped_reads/Progne_murphyi_Peru_LSUMZ_114186.fxm.sorted.md.bam 
    Progne_sinaloae_Mexico_KU_40044: /data3/hirundinidae_phylogeny/workflow_reanalysis/phyluce/1-mapping/mapped_reads/Progne_sinaloae_Mexico_KU_40044.fxm.sorted.md.bam 
    Progne_sinaloae_Mexico_KU_40045: /data3/hirundinidae_phylogeny/workflow_reanalysis/phyluce/1-mapping/mapped_reads/Progne_sinaloae_Mexico_KU_40045.fxm.sorted.md.bam 
    Psalidoprocne_obscura_Ivory_Coast_FMNH_279082: /data3/hirundinidae_phylogeny/workflow_reanalysis/phyluce/1-mapping/mapped_reads/Psalidoprocne_obscura_Ivory_Coast_FMNH_279082.fxm.sorted.md.bam 
    Pseudochelidon_eurystomina_Congo_FMNH_213526: /data3/hirundinidae_phylogeny/workflow_reanalysis/phyluce/1-mapping/mapped_reads/Pseudochelidon_eurystomina_Congo_FMNH_213526.fxm.sorted.md.bam 
    Pseudochelidon_sirintarae_Thailand_AMNH_SKIN_708673: /data3/hirundinidae_phylogeny/workflow_reanalysis/phyluce/1-mapping/mapped_reads/Pseudochelidon_sirintarae_Thailand_AMNH_SKIN_708673.fxm.sorted.md.bam 
    Ptyonoprogne_concolor_India_FMNH_233366: /data3/hirundinidae_phylogeny/workflow_reanalysis/phyluce/1-mapping/mapped_reads/Ptyonoprogne_concolor_India_FMNH_233366.fxm.sorted.md.bam 
    Ptyonoprogne_obsoleta_Egypt_UMMZ_224076: /data3/hirundinidae_phylogeny/workflow_reanalysis/phyluce/1-mapping/mapped_reads/Ptyonoprogne_obsoleta_Egypt_UMMZ_224076.fxm.sorted.md.bam 
    Riparia_congica_Congo_FMNH_213543: /data3/hirundinidae_phylogeny/workflow_reanalysis/phyluce/1-mapping/mapped_reads/Riparia_congica_Congo_FMNH_213543.fxm.sorted.md.bam 
contigs:
    Cecropis_daurica_Africa_Uganda_LACM_71522: /data3/hirundinidae_phylogeny/Hirundinidae_UCE/phyluce/spades-assemblies-ingroup-skins/Cecropis_daurica_Africa_Uganda_LACM_71522_spades/contigs.fasta
    Delichon_nipalense_Nepal_FMNH_276821: /data3/hirundinidae_phylogeny/Hirundinidae_UCE/phyluce/spades-assemblies-ingroup-skins/Delichon_nipalense_Nepal_FMNH_276821_spades/contigs.fasta
    Hirundo_lucida_Gambia_FMNH_320046: /data3/hirundinidae_phylogeny/Hirundinidae_UCE/phyluce/spades-assemblies-ingroup-skins/Hirundo_lucida_Gambia_FMNH_320046_spades/contigs.fasta
    Hirundo_megaensis_Ethiopia_AMNH_SKIN_348889: /data3/hirundinidae_phylogeny/Hirundinidae_UCE/phyluce/spades-assemblies-ingroup-skins/Hirundo_megaensis_Ethiopia_AMNH_SKIN_348889_spades/contigs.fasta
    Hirundo_nigrorufa_Angola_AMNH_SKIN_707943: /data3/hirundinidae_phylogeny/Hirundinidae_UCE/phyluce/spades-assemblies-ingroup-skins/Hirundo_nigrorufa_Angola_AMNH_SKIN_707943_spades/contigs.fasta
    Petrochelidon_fluvicola_India_FMNH_233414: /data3/hirundinidae_phylogeny/Hirundinidae_UCE/phyluce/spades-assemblies-ingroup-skins/Petrochelidon_fluvicola_India_FMNH_233414_spades/contigs.fasta
    Petrochelidon_rufigula_Angola_AMNH_SKIN_707945: /data3/hirundinidae_phylogeny/Hirundinidae_UCE/phyluce/spades-assemblies-ingroup-skins/Petrochelidon_rufigula_Angola_AMNH_SKIN_707945_spades/contigs.fasta
    Phedinopsis_brazzae_Angola_AMNH_SKIN_764767: /data3/hirundinidae_phylogeny/Hirundinidae_UCE/phyluce/spades-assemblies-ingroup-skins/Phedinopsis_brazzae_Angola_AMNH_SKIN_764767_spades/contigs.fasta
    Progne_cryptoleuca_Cuba_LSUMZ_142940: /data3/hirundinidae_phylogeny/Hirundinidae_UCE/phyluce/spades-assemblies-ingroup-skins/Progne_cryptoleuca_Cuba_LSUMZ_142940_spades/contigs.fasta
    Progne_cryptoleuca_Cuba_LSUMZ_142944: /data3/hirundinidae_phylogeny/Hirundinidae_UCE/phyluce/spades-assemblies-ingroup-skins/Progne_cryptoleuca_Cuba_LSUMZ_142944_spades/contigs.fasta
    Progne_modesta_Ecuador_MVZ_130129: /data3/hirundinidae_phylogeny/Hirundinidae_UCE/phyluce/spades-assemblies-ingroup-skins/Progne_modesta_Ecuador_MVZ_130129_spades/contigs.fasta
    Progne_murphyi_Peru_LSUMZ_114185: /data3/hirundinidae_phylogeny/Hirundinidae_UCE/phyluce/spades-assemblies-ingroup-skins/Progne_murphyi_Peru_LSUMZ_114185_spades/contigs.fasta
    Progne_murphyi_Peru_LSUMZ_114186: /data3/hirundinidae_phylogeny/Hirundinidae_UCE/phyluce/spades-assemblies-ingroup-skins/Progne_murphyi_Peru_LSUMZ_114186_spades/contigs.fasta
    Progne_sinaloae_Mexico_KU_40044: /data3/hirundinidae_phylogeny/Hirundinidae_UCE/phyluce/spades-assemblies-ingroup-skins/Progne_sinaloae_Mexico_KU_40044_spades/contigs.fasta
    Progne_sinaloae_Mexico_KU_40045: /data3/hirundinidae_phylogeny/Hirundinidae_UCE/phyluce/spades-assemblies-ingroup-skins/Progne_sinaloae_Mexico_KU_40045_spades/contigs.fasta
    Psalidoprocne_obscura_Ivory_Coast_FMNH_279082: /data3/hirundinidae_phylogeny/Hirundinidae_UCE/phyluce/spades-assemblies-ingroup-skins/Psalidoprocne_obscura_Ivory_Coast_FMNH_279082_spades/contigs.fasta
    Pseudochelidon_eurystomina_Congo_FMNH_213526: /data3/hirundinidae_phylogeny/Hirundinidae_UCE/phyluce/spades-assemblies-ingroup-skins/Pseudochelidon_eurystomina_Congo_FMNH_213526_spades/contigs.fasta
    Pseudochelidon_sirintarae_Thailand_AMNH_SKIN_708673: /data3/hirundinidae_phylogeny/Hirundinidae_UCE/phyluce/spades-assemblies-ingroup-skins/Pseudochelidon_sirintarae_Thailand_AMNH_SKIN_708673_spades/contigs.fasta
    Ptyonoprogne_concolor_India_FMNH_233366: /data3/hirundinidae_phylogeny/Hirundinidae_UCE/phyluce/spades-assemblies-ingroup-skins/Ptyonoprogne_concolor_India_FMNH_233366_spades/contigs.fasta
    Ptyonoprogne_obsoleta_Egypt_UMMZ_224076: /data3/hirundinidae_phylogeny/Hirundinidae_UCE/phyluce/spades-assemblies-ingroup-skins/Ptyonoprogne_obsoleta_Egypt_UMMZ_224076_spades/contigs.fasta
    Riparia_congica_Congo_FMNH_213543: /data3/hirundinidae_phylogeny/Hirundinidae_UCE/phyluce/spades-assemblies-ingroup-skins/Riparia_congica_Congo_FMNH_213543_spades/contigs.fasta
```

#### Run correction workflow

```
phyluce_workflow --config ./2-correction.config --output 2-correction --workflow correction --cores 24
```

#### Output files

The outputs from the correction step are in:
```
./2-correction/consensus
./2-correction/filtered_norm_pileups
```

## Part 2 - Finding UCE loci in the corrected skin sample contigs

We now have a set of corrected consensus contig .fasta sequences for the 21 skin samples.

We need to identify UCE loci among the contigs, to be integrated into datasets with the tissue samples.

------------------------------------------------------------------------------------------
### Input data locations

The corrected contigs from the correction step are in:
```
/data3/hirundinidae_phylogeny/workflow_reanalysis/phyluce/2-correction/consensus
```

The 5k UCE probe set is in:
```
/data3/hirundinidae_phylogeny/Hirundinidae_UCE/phyluce/uce-5k-probes.fasta
```

The probe set file is small. For convenience/posterity, retrieve it to the working directory.
```
cd /data3/hirundinidae_phylogeny/workflow_reanalysis/phyluce
cp /data3/hirundinidae_phylogeny/Hirundinidae_UCE/phyluce/uce-5k-probes.fasta .
```

------------------------------------------------------------------------------------------
### 1. Format corrected consensus contig file names

We want these to match the inputs for the outgroup/tissue samples.

Currently, the filename extensions are `<sample>.consensus.filt.fasta`. We want them to be `<sample>.contigs.fasta`.

We'll write copies to a separate `consensus-rename` directory, to which we'll also add the .fasta files for the outgroup/tissue samples in the next step.

```
cd ./2-correction/
mkdir consensus-rename
cp ./consensus/*.fasta ./consensus-rename/
cd consensus-rename
for i in *.fasta; do name=`echo $i | cut -d'.' -f1`; echo $name; mv ./$name.consensus.filt.fasta ./$name.contigs.fasta; done
```

------------------------------------------------------------------------------------------
### 2. Retrieve contigs for outgroup/tissue samples

We will want to find UCE loci for all of the samples together so that they have a common database.
```
cd /data3/hirundinidae_phylogeny/workflow_reanalysis/phyluce/2-correction/consensus-rename
cp /data3/hirundinidae_phylogeny/Hirundinidae_UCE/phyluce/spades-assemblies-outgroup/contigs/*.fasta .
cp /data3/hirundinidae_phylogeny/Hirundinidae_UCE/phyluce/spades-assemblies-ingroup-tissues/contigs/*.fasta .
cp /data3/hirundinidae_phylogeny/Hirundinidae_UCE/phyluce/spades-assemblies-tachycineta-project/contigs/*.fasta .
```

We now have the data for the total 122 samples!

------------------------------------------------------------------------------------------
### 3. Find UCE loci in corrected consensus contigs + assembled contigs for outgroups/ingroup tissue samples.

Run `phyluce_assembly_match_contigs_to_probes` to find UCE matches.
```
conda activate phyluce
phyluce_assembly_match_contigs_to_probes --contigs 2-correction/consensus-rename/ --probes uce-5k-probes.fasta	--output 3-uce-search-results
```

------------------------------------------------------------------------------------------
### 4. Summary of results

The output probe matching stats for the skin samples are in `/data3/hirundinidae_phylogeny/workflow_reanalysis/Hirundinidae_UCE_skins_reanalysis.xlsx`.

Outputs are in:
```
./3-uce-search-results/<sample>.contigs.lastz
./3-uce-search-results/probe.matches.sqlite
```

## Part 3 - Extract UCE loci for skins; generate datasets with tissue samples

We have matched the corrected consensus contigs for the 21 skin samples to the 5k UCE probes, together with sequences for the outgroups and ingroup tissue samples.

We'll now extract UCE sequences and generate locus datasets across the combined set of samples.

We'll run several phyluce commands in this part of the workflow to extract/format UCE loci.

------------------------------------------------------------------------------------------
### Input data locations

The contigs uniquely matched to UCE probes are in:
```
/data3/hirundinidae_phylogeny/workflow_reanalysis/phyluce/3-uce-search-results
```

------------------------------------------------------------------------------------------
### Set up environment

```
cd /data3/hirundinidae_phylogeny/workflow_reanalysis/phyluce/
mkdir 4-datasets
cd 4-datasets
mkdir all
```

------------------------------------------------------------------------------------------
### 1. Run `phyluce_assembly_get_match_counts` on combined dataset

#### Format dataset configuration file

```
/data3/hirundinidae_phylogeny/workflow_reanalysis/phyluce/4-datasets.config
```

Note: some (13) of the samples have repeated Country, Locality names (e.g., Cecropis_daurica_China_China_KU_99748). Clare addressed this in her '7d-rename assemblies' document, but I'll proceed with the names as is.

#### Format taxon-locus configuration file

```
/data3/hirundinidae_phylogeny/workflow_reanalysis/phyluce/4-datasets/all/all-samples-incomplete.config

```

#### Run analysis
```
cd /data3/hirundinidae_phylogeny/workflow_reanalysis/phyluce/
conda activate phyluce
phyluce_assembly_get_match_counts --locus-db ./3-uce-search-results/probe.matches.sqlite --taxon-list-config ./4-datasets.config --taxon-group 'all' --incomplete-matrix --output ./4-datasets/all/all-samples-incomplete.config --log-path .
```

------------------------------------------------------------------------------------------
### 2. Run `phyluce_assembly_get_fastas_from_match_counts` to format .fasta files from the matched UCE contigs

This command references outputs from the last few steps, including the assembled/corrected contigs, search matches, etc.

```
phyluce_assembly_get_fastas_from_match_counts --contigs ./2-correction/consensus-rename/ --locus-db ./3-uce-search-results/probe.matches.sqlite --match-count-output ./4-datasets/all/all-samples-incomplete.config --output ./4-datasets/all/all-samples-incomplete.fasta --incomplete-matrix ./4-datasets/all/all-samples-incomplete.incomplete --log-path .
```

------------------------------------------------------------------------------------------
### 3. Run `phyluce_assembly_explode_get_fastas_file` to generate taxon-specific fasta files

This will produce a single fasta file for each taxon, including all of the UCEs for the sample.

```
phyluce_assembly_explode_get_fastas_file --input ./4-datasets/all/all-samples-incomplete.fasta --output ./4-datasets/all/exploded-fastas-all-samples-incomplete --by-taxon
```

------------------------------------------------------------------------------------------
### 4. Run `phyluce_assembly_get_fasta_lengths` to generate summary statistics on extracted UCE loci

```
for seq in ./4-datasets/all/exploded-fastas-all-samples-incomplete/*.fasta; do phyluce_assembly_get_fasta_lengths --input $seq --csv; done > ./4-datasets/all/all-samples-incomplete-stats.txt
```

## Part 4 - Align UCE loci

We've extracted UCE sequences for the corrected skin samples and the outgroup/tissue samples in the previous steps.

Now we'll align and format the data for analysis.

------------------------------------------------------------------------------------------
### Input data

Fasta sequences for the all-taxon dataset are in `/data3/hirundinidae_phylogeny/workflow_reanalysis/phyluce/4-datasets/all/`

------------------------------------------------------------------------------------------
### 1. Use phyluce to run MAFFT alignment on input data

Note: we'll specify `--ambiguous` to not remove sequences with ambiguous bases (i.e., introduced in the correction step for skin samples).

```
cd /data3/hirundinidae_phylogeny/workflow_reanalysis/phyluce/
conda activate phyluce
phyluce_align_seqcap_align --input ./4-datasets/all/all-samples-incomplete.fasta --output ./4-datasets/all/mafft-fasta-untrimmed-all-samples-incomplete --taxa 122 --aligner mafft --cores 24 --incomplete-matrix --output-format fasta --no-trim --ambiguous --log-path .
```

Note: loci uce-2821, uce-3732, uce-4427, uce-7305, uce-2241, uce-5887, uce-4726, uce-5673 dropped due to too few taxa (n < 3).

------------------------------------------------------------------------------------------
### 2. Rename sequences in alignments

We need the sequence headers to be compatible with Gblocks, so will shorten to just the museum identifier.
```
cd /data3/hirundinidae_phylogeny/workflow_reanalysis/phyluce/4-datasets/all
cp -r mafft-fasta-untrimmed-all-samples-incomplete/ mafft-fasta-untrimmed-all-samples-incomplete-SHORTNAMES
```

Format original and matched short name files.
```
all_originalNames.txt
all_shortNames.txt
```

Note: the original names differ slightly for some samples than Clare's file, because the small set of samples have the double locality, country parts of the name. The details of the conversion are in the 'Hirundinidae_taxon_name_conversion.xlsx' file.

Run the conversion, then remove the temporary fasta files.
```
for file in mafft-fasta-untrimmed-all-samples-incomplete-SHORTNAMES/*.fasta; do echo $file >> ./rename_fasta_all_samples_incomplete.log; for i in {1..122}; do sed -i.temp.fasta "s/$(head -$i all_originalNames.txt | tail -1)/$(head -$i all_shortNames.txt | tail -1)/" $file; done; done
rm ./mafft-fasta-untrimmed-all-samples-incomplete-SHORTNAMES/*.temp.fasta

```

------------------------------------------------------------------------------------------
### 3. Run `phyluce_align_get_gblocks_trimmed_alignments_from_untrimmed` to trim alignments

This calls Gblocks to trim alignments based on several input parameters.
```
cd /data3/hirundinidae_phylogeny/workflow_reanalysis/phyluce/
mkdir temp_log
phyluce_align_get_gblocks_trimmed_alignments_from_untrimmed --alignments ./4-datasets/all/mafft-fasta-untrimmed-all-samples-incomplete-SHORTNAMES --output ./4-datasets/all/mafft-gblocks-nexus-internal-trimmed-all-samples-incomplete-SHORTNAMES --b2 0.65 --b4 8 --input-format fasta --output-format nexus --cores 24 --log ./temp_log
```

------------------------------------------------------------------------------------------
### 4. Rename sequences in aligned/trimmed nexus files with original names

```
cd ./4-datasets/all/
cp -r mafft-gblocks-nexus-internal-trimmed-all-samples-incomplete-SHORTNAMES mafft-gblocks-nexus-internal-trimmed-all-samples-incomplete

for file in mafft-gblocks-nexus-internal-trimmed-all-samples-incomplete/*.nexus;
	do echo $file >> ./rename_nexus_all_samples_incomplete.log;
	for i in {1..122};
		do sed -i.temp.nexus "s/$(head -$i all_shortNames.txt | tail -1)/$(head -$i all_originalNames.txt | tail -1)/" $file;
	done;
done

rm ./mafft-gblocks-nexus-internal-trimmed-all-samples-incomplete/*.temp.nexus
```

------------------------------------------------------------------------------------------
### 5. Fix issue with weird double name headers

For some reason, the steps above produced a subset of nexus files with doubled names.

E.g., 'Alopochelidon_fucata_Uruguay_CUMV_50652' and 'Alopochelidon_fucata_Uruguay_Alopochelidon_fucata_Uruguay_CUMV_50652'.

This presents a big problem later on, when concatenating the sequences together, as the two different versions are seen as two separate taxa.

### Format list of double names to search nexus files with

```
all_doubleNames.txt
```

#### Do grep search of all nexus files in `./mafft-gblocks-nexus-internal-trimmed-all-samples-incomplete/`

Do search and output to temporary file:
```
for i in `cat all_doubleNames.txt`; do for nex in ./mafft-gblocks-nexus-internal-trimmed-all-samples-incomplete/*.nexus; do grep $i $nex; done | cut -d"'" -f2 | cut -d'_' -f1 | uniq >> all_doubleNames.nexus.tmp.txt; done
```

Get unique nexus file hits:
```
awk '{print $1}' all_doubleNames.nexus.tmp.txt | sort | uniq > all_doubleNames.nexus.txt
```

#### Run sed replace script on nexus files with messed up sample names

Wrote `fixNexusDoubleNames.sh` to call list above and replace with correct names.
```
sh fixNexusDoubleNames.sh
```

Note: this relies on the input list being `./all_doubleNames.nexus.txt` and the target directory `./mafft-gblocks-nexus-internal-trimmed-all-samples-incomplete/`.

#### Check that the script worked

There should not be any of the double names in any of the nexus files.

```
for nex in ./mafft-gblocks-nexus-internal-trimmed-all-samples-incomplete/*.nexus; do grep 'Alopochelidon_fucata_Uruguay_Alopochelidon_fucata_Uruguay_CUMV_50652' $nex; done | cut -d"'" -f2 | cut -d'_' -f1 | uniq
```

This seems to have done the job.

------------------------------------------------------------------------------------------
### 6. Get alignment summary statistics

Run `phyluce_align_get_align_summary_data` to summarize stats on informative sites, lengths, missingness, etc.
```
cd /data3/hirundinidae_phylogeny/workflow_reanalysis/phyluce/
phyluce_align_get_align_summary_data --alignments ./4-datasets/all/mafft-gblocks-nexus-internal-trimmed-all-samples-incomplete --cores 24 --log-path ./temp_log
```

------------------------------------------------------------------------------------------
### 7. Remove locus names from nexus file header lines

```
phyluce_align_remove_locus_name_from_files --alignments ./4-datasets/all/mafft-gblocks-nexus-internal-trimmed-all-samples-incomplete --output ./4-datasets/all/mafft-gblocks-nexus-internal-trimmed-all-samples-incomplete-clean --cores 24 --log-path ./temp_log
```

## Part 5 - Format data matrices for analysis

We've aligned UCE loci using mafft and Gblocks and formatted nexus alignments for analysis.

Now we need to format input datasets for partitioning and analysis in RAxML, SVDquartets, Astral, etc.

------------------------------------------------------------------------------------------
### 1. Format dataset for all samples for pilot analyses

We'll use the results from this dataset to determine if any samples with low quality input data need to be removed from downstream analysis.

#### Format UCE locus dataset with 95% samples per locus

Run `phyluce_align_get_only_loci_with_min_taxa` to extract 95% complete matrix.
```
cd /data3/hirundinidae_phylogeny/workflow_renalysis/phyluce
conda activate phyluce
phyluce_align_get_only_loci_with_min_taxa --alignments ./4-datasets/all/mafft-gblocks-nexus-internal-trimmed-all-samples-incomplete-clean --output ./4-datasets/all/mafft-gblocks-nexus-internal-trimmed-all-samples-incomplete-clean-95p --taxa 122 --percent 0.95 --cores 24 --log-path ./temp_log
```

This retains 3442 of 5026 total alignments with >= 95% of taxa (n = 122 taxa).

------------------------------------------------------------------------------------------
### 2. Format revised dataset 1, with highly problematic skin samples removed

Following preliminary analysis, it was clear that 4 samples should be removed from analysis (see 'README-7a-RAxML-concatenated-analysis-pilot-notes.docx').

Here, we'll remove those samples from the alignments and generate a 95% complete data matrix for downstream analysis.

We can effectively consider this the 'full' dataset moving forward.

#### Set up environment

```
cd /data3/hirundinidae_phylogeny/workflow_renalysis/phyluce/4-datasets/
mkdir all-rev1
cd ..
```

#### Extract alignments without excluded samples

Extract alignments:
```
phyluce_align_extract_taxa_from_alignments --alignments ./4-datasets/all/mafft-gblocks-nexus-internal-trimmed-all-samples-incomplete-clean --output ./4-datasets/all-rev1/mafft-gblocks-nexus-internal-trimmed-all-rev1-samples-incomplete-clean --output-format nexus --exclude Progne_cryptoleuca_Cuba_LSUMZ_142940 Progne_cryptoleuca_Cuba_LSUMZ_142944 Riparia_congica_Congo_FMNH_213543 Progne_modesta_Ecuador_MVZ_130129 --cores 8 --log-path ./temp_log
```

Verify that the excluded taxa are not in the resulting alignments:
```
for i in ./4-datasets/all-rev1/mafft-gblocks-nexus-internal-trimmed-all-rev1-samples-incomplete-clean/*.nexus; do grep 'Progne_cryptoleuca_Cuba_LSUMZ_142940' $i; done
for i in ./4-datasets/all-rev1/mafft-gblocks-nexus-internal-trimmed-all-rev1-samples-incomplete-clean/*.nexus; do grep 'Progne_cryptoleuca_Cuba_LSUMZ_142944' $i; done
for i in ./4-datasets/all-rev1/mafft-gblocks-nexus-internal-trimmed-all-rev1-samples-incomplete-clean/*.nexus; do grep 'Riparia_congica_Congo_FMNH_213543' $i; done
for i in ./4-datasets/all-rev1/mafft-gblocks-nexus-internal-trimmed-all-rev1-samples-incomplete-clean/*.nexus; do grep 'Progne_modesta_Ecuador_MVZ_130129' $i; done
```

#### Format UCE locus dataset with 95% samples per locus

```
phyluce_align_get_only_loci_with_min_taxa --alignments ./4-datasets/all-rev1/mafft-gblocks-nexus-internal-trimmed-all-rev1-samples-incomplete-clean --output ./4-datasets/all-rev1/mafft-gblocks-nexus-internal-trimmed-all-rev1-samples-incomplete-clean-95p --taxa 118 --percent 0.95 --cores 8 --log-path ./temp_log
```

This retains 4009 of 5026 total alignments with 95% of taxa (n = 118 taxa).

------------------------------------------------------------------------------------------
### 3. Format revised dataset 2, with only tissue samples + Pseudochelidon samples

Here, we'll remove skin samples from the alignments, generate a 95% complete data matrix, and proceed with analyses focusing on only these samples.

#### Set up environment

```
cd /data3/hirundinidae_phylogeny/workflow_renalysis/phyluce/4-datasets/
mkdir all-rev2
cd ..
```

#### Extract alignments without excluded samples

Extract alignments:
```
phyluce_align_extract_taxa_from_alignments --alignments ./4-datasets/all/mafft-gblocks-nexus-internal-trimmed-all-samples-incomplete-clean --output ./4-datasets/all-rev2/mafft-gblocks-nexus-internal-trimmed-all-rev2-samples-incomplete-clean --output-format nexus --exclude Cecropis_daurica_Africa_Uganda_LACM_71522 Delichon_nipalense_Nepal_FMNH_276821 Hirundo_lucida_Gambia_FMNH_320046 Hirundo_megaensis_Ethiopia_AMNH_SKIN_348889 Hirundo_nigrorufa_Angola_AMNH_SKIN_707943 Petrochelidon_fluvicola_India_FMNH_233414 Petrochelidon_rufigula_Angola_AMNH_SKIN_707945 Phedinopsis_brazzae_Angola_AMNH_SKIN_764767 Progne_cryptoleuca_Cuba_LSUMZ_142940 Progne_cryptoleuca_Cuba_LSUMZ_142944 Progne_modesta_Ecuador_MVZ_130129 Progne_murphyi_Peru_LSUMZ_114185 Progne_murphyi_Peru_LSUMZ_114186 Progne_sinaloae_Mexico_KU_40044 Progne_sinaloae_Mexico_KU_40045 Psalidoprocne_obscura_Ivory_Coast_FMNH_279082 Ptyonoprogne_concolor_India_FMNH_233366 Ptyonoprogne_obsoleta_Egypt_UMMZ_224076 Riparia_congica_Congo_FMNH_213543 --cores 8 --log-path ./temp_log
```

Verify that the excluded taxa are not in the resulting alignments:
```
for i in ./4-datasets/all-rev2/mafft-gblocks-nexus-internal-trimmed-all-rev2-samples-incomplete-clean/*.nexus; do grep 'Progne_cryptoleuca_Cuba_LSUMZ_142940' $i; done
for i in ./4-datasets/all-rev2/mafft-gblocks-nexus-internal-trimmed-all-rev2-samples-incomplete-clean/*.nexus; do grep 'Cecropis_daurica_Africa_Uganda_LACM_71522' $i; done
for i in ./4-datasets/all-rev2/mafft-gblocks-nexus-internal-trimmed-all-rev2-samples-incomplete-clean/*.nexus; do grep 'Phedinopsis_brazzae_Angola_AMNH_SKIN_764767' $i; done
for i in ./4-datasets/all-rev2/mafft-gblocks-nexus-internal-trimmed-all-rev2-samples-incomplete-clean/*.nexus; do grep 'Progne_modesta_Ecuador_MVZ_130129' $i; done
```

#### Format UCE locus dataset with 95% samples per locus

```
phyluce_align_get_only_loci_with_min_taxa --alignments ./4-datasets/all-rev2/mafft-gblocks-nexus-internal-trimmed-all-rev2-samples-incomplete-clean --output ./4-datasets/all-rev2/mafft-gblocks-nexus-internal-trimmed-all-rev2-samples-incomplete-clean-95p --taxa 103 --percent 0.95 --cores 8 --log-path ./temp_log
```

This retains 4518 of 5001 total alignments with 95% of taxa (n = 103 taxa).

------------------------------------------------------------------------------------------
### 4. Format revised dataset 3, with only tissue samples

Here, we'll remove all skin samples from the alignments, generate a 95% complete data matrix, and proceed with analyses focusing on only these samples.

#### Set up environment

```
cd /data3/hirundinidae_phylogeny/workflow_renalysis/phyluce/4-datasets/
mkdir all-rev3
cd ..
```

#### Extract alignments without excluded samples

Extract alignments:
```
phyluce_align_extract_taxa_from_alignments --alignments ./4-datasets/all/mafft-gblocks-nexus-internal-trimmed-all-samples-incomplete-clean --output ./4-datasets/all-rev3/mafft-gblocks-nexus-internal-trimmed-all-rev3-samples-incomplete-clean --output-format nexus --exclude Pseudochelidon_eurystomina_Congo_FMNH_213526 Pseudochelidon_sirintarae_Thailand_AMNH_SKIN_708673 Cecropis_daurica_Africa_Uganda_LACM_71522 Delichon_nipalense_Nepal_FMNH_276821 Hirundo_lucida_Gambia_FMNH_320046 Hirundo_megaensis_Ethiopia_AMNH_SKIN_348889 Hirundo_nigrorufa_Angola_AMNH_SKIN_707943 Petrochelidon_fluvicola_India_FMNH_233414 Petrochelidon_rufigula_Angola_AMNH_SKIN_707945 Phedinopsis_brazzae_Angola_AMNH_SKIN_764767 Progne_cryptoleuca_Cuba_LSUMZ_142940 Progne_cryptoleuca_Cuba_LSUMZ_142944 Progne_modesta_Ecuador_MVZ_130129 Progne_murphyi_Peru_LSUMZ_114185 Progne_murphyi_Peru_LSUMZ_114186 Progne_sinaloae_Mexico_KU_40044 Progne_sinaloae_Mexico_KU_40045 Psalidoprocne_obscura_Ivory_Coast_FMNH_279082 Ptyonoprogne_concolor_India_FMNH_233366 Ptyonoprogne_obsoleta_Egypt_UMMZ_224076 Riparia_congica_Congo_FMNH_213543 --cores 8 --log-path ./temp_log
```

Verify that the excluded taxa are not in the resulting alignments:
```
for i in ./4-datasets/all-rev2/mafft-gblocks-nexus-internal-trimmed-all-rev2-samples-incomplete-clean/*.nexus; do grep 'Progne_cryptoleuca_Cuba_LSUMZ_142940' $i; done
for i in ./4-datasets/all-rev2/mafft-gblocks-nexus-internal-trimmed-all-rev2-samples-incomplete-clean/*.nexus; do grep 'Cecropis_daurica_Africa_Uganda_LACM_71522' $i; done
for i in ./4-datasets/all-rev2/mafft-gblocks-nexus-internal-trimmed-all-rev2-samples-incomplete-clean/*.nexus; do grep 'Phedinopsis_brazzae_Angola_AMNH_SKIN_764767' $i; done
for i in ./4-datasets/all-rev2/mafft-gblocks-nexus-internal-trimmed-all-rev2-samples-incomplete-clean/*.nexus; do grep 'Progne_modesta_Ecuador_MVZ_130129' $i; done
```

#### Format UCE locus dataset with 95% samples per locus

```
phyluce_align_get_only_loci_with_min_taxa --alignments ./4-datasets/all-rev3/mafft-gblocks-nexus-internal-trimmed-all-rev3-samples-incomplete-clean --output ./4-datasets/all-rev3/mafft-gblocks-nexus-internal-trimmed-all-rev3-samples-incomplete-clean-95p --taxa 101 --percent 0.95 --cores 8 --log-path ./temp_log
```

This retains 4565 of 4986 total alignments with 95% of taxa (n = 101 taxa).

Get alignment summary data:
```
phyluce_align_get_align_summary_data --alignments ./4-datasets/all-rev3/mafft-gblocks-nexus-internal-trimmed-all-rev3-samples-incomplete-clean-95p --cores 24 --log-path ./temp_log
```

## Part 6 - Concatenate UCE alignments

We now have a 95% complete UCE data matrix for all samples (i.e., outgroups, tissue, and skin samples).

We will run a concatenated maximum likelihood analysis in RAxML and a coalescent-based analysis in SVDquartets. To do this we need a concatenated input alignment.

### Set up environment

```
cd /data3/hirundinidae_phylogeny/workflow_renalysis/phyluce
mkdir 5-input-alignments
cd 5-input-alignments
mkdir all
mkdir all-rev1
``` 

------------------------------------------------------------------------------------------
### 1. Full dataset (all 122 taxa)

#### 1. Use phyluce to generate concatenated alignment nexus

```
cd /data3/hirundinidae_phylogeny/workflow_renalysis/phyluce
conda activate phyluce
phyluce_align_concatenate_alignments --alignments ./4-datasets/all/mafft-gblocks-nexus-internal-trimmed-all-samples-incomplete-clean-95p --output ./5-input-alignments/all/hirundinidae-all-95-nexus-concat --nexus --log-path .
```

Confirm that there are 122 taxon labels in the nexus:
```
cut -d' ' -f1 ./5-input-alignments/all/hirundinidae-all-95-nexus-concat/hirundinidae-all-95-nexus-concat.nexus
```

Yes, there are.

#### 2. Use phyluce to generate concatenated alignment phylip

```
phyluce_align_concatenate_alignments --alignments ./4-datasets/all/mafft-gblocks-nexus-internal-trimmed-all-samples-incomplete-clean-95p --output ./5-input-alignments/all/hirundinidae-all-95-phylip-concat --phylip --log-path .
```
The concatenated phylip alignment is in `./5-input-alignments/all/hirundinidae-all-95-phylip-concat/hirundinidae-all-95-phylip-concat.phylip`.

#### 3. Format UCE labels in concatenated outputs

We want 'uce_XXXX' instead of 'uce-XXXX' format.

```
sed -i 's/uce-/uce_/g' ./5-input-alignments/all/hirundinidae-all-95-nexus-concat/hirundinidae-all-95-nexus-concat.nexus
sed -i 's/uce-/uce_/g' ./5-input-alignments/all/hirundinidae-all-95-phylip-concat/hirundinidae-all-95-phylip-concat.charsets
```

#### 4. Format concatenated alignment for pilot dataset

To diagnose how 'successful' the phyluce correction step was for the skin samples, it'll be good to run a pilot analysis on a subset of loci to determine if the long branch length issue was solved.

##### 1. Select random 100 UCE loci for analysis

```
cd /data3/hirundinidae_phylogeny/workflow_renalysis/phyluce
mkdir ./4-datasets/all/mafft-gblocks-nexus-internal-trimmed-all-samples-incomplete-clean-95p-pilot
cd ./4-datasets/all/mafft-gblocks-nexus-internal-trimmed-all-samples-incomplete-clean-95p
list=`ls -lh *.nexus | cut -d':' -f2 | cut -d' ' -f2 | shuf -n 100`; for i in $list; do cp $i ../mafft-gblocks-nexus-internal-trimmed-all-samples-incomplete-clean-95p-pilot; done
```

##### 2. Use phyluce to generate concatenated alignment phylip

```
cd /data3/hirundinidae_phylogeny/workflow_renalysis/phyluce
conda activate phyluce
phyluce_align_concatenate_alignments --alignments ./4-datasets/all/mafft-gblocks-nexus-internal-trimmed-all-samples-incomplete-clean-95p-pilot --output ./5-input-alignments/all/hirundinidae-all-95-pilot-phylip-concat --phylip --log-path .
```

##### 3. Format UCE labels

```
sed -i 's/uce-/uce_/g' ./5-input-alignments/all/hirundinidae-all-95-pilot-phylip-concat/hirundinidae-all-95-pilot-phylip-concat.charsets
```

------------------------------------------------------------------------------------------
### 2. Full, revised 1 dataset (n = 118 taxa)

Here, only the 4 most problematic skin samples have been removed.

#### 1. Format concatenated alignment for pilot dataset (100 loci, as above)

##### 1. Select random 100 UCE loci for analysis

```
cd /data3/hirundinidae_phylogeny/workflow_reanalysis/phyluce
mkdir ./4-datasets/all-rev1/mafft-gblocks-nexus-internal-trimmed-all-rev1-samples-incomplete-clean-95p-pilot
cd ./4-datasets/all-rev1/mafft-gblocks-nexus-internal-trimmed-all-rev1-samples-incomplete-clean-95p
list=`ls -lh *.nexus | cut -d':' -f2 | cut -d' ' -f2 | shuf -n 100`; for i in $list; do cp $i ../mafft-gblocks-nexus-internal-trimmed-all-rev1-samples-incomplete-clean-95p-pilot; done
```

##### 2. Use phyluce to generate concatenated alignment phylip

```
cd /data3/hirundinidae_phylogeny/workflow_reanalysis/phyluce
phyluce_align_concatenate_alignments --alignments ./4-datasets/all-rev1/mafft-gblocks-nexus-internal-trimmed-all-rev1-samples-incomplete-clean-95p-pilot --output ./5-input-alignments/all-rev1/hirundinidae-all-rev1-95-pilot-phylip-concat --phylip --log-path .
```

##### 3. Use phyluce to generate concatenated alignment nexus

```
phyluce_align_concatenate_alignments --alignments ./4-datasets/all-rev1/mafft-gblocks-nexus-internal-trimmed-all-rev1-samples-incomplete-clean-95p-pilot --output ./5-input-alignments/all-rev1/hirundinidae-all-rev1-95-pilot-nexus-concat --nexus --log-path .
```

##### 4. Format UCE labels

```
sed -i 's/uce-/uce_/g' ./5-input-alignments/all-rev1/hirundinidae-all-rev1-95-pilot-phylip-concat/hirundinidae-all-rev1-95-pilot-phylip-concat.charsets
sed -i 's/uce-/uce_/g' ./5-input-alignments/all-rev1/hirundinidae-all-rev1-95-pilot-nexus-concat/hirundinidae-all-rev1-95-pilot-nexus-concat.nexus
```

#### 2. Format concatenated alignment for full dataset

##### 1. Use phyluce to generate concatenated alignment phylip

```
cd /data3/hirundinidae_phylogeny/workflow_reanalysis/phyluce
phyluce_align_concatenate_alignments --alignments ./4-datasets/all-rev1/mafft-gblocks-nexus-internal-trimmed-all-rev1-samples-incomplete-clean-95p --output ./5-input-alignments/all-rev1/hirundinidae-all-rev1-95-phylip-concat --phylip --log-path .
```

##### 2. Use phyluce to generate concatenated alignment nexus

```
phyluce_align_concatenate_alignments --alignments ./4-datasets/all-rev1/mafft-gblocks-nexus-internal-trimmed-all-rev1-samples-incomplete-clean-95p --output ./5-input-alignments/all-rev1/hirundinidae-all-rev1-95-nexus-concat --nexus --log-path .
```

##### 3. Format UCE labels

```
sed -i 's/uce-/uce_/g' ./5-input-alignments/all-rev1/hirundinidae-all-rev1-95-phylip-concat/hirundinidae-all-rev1-95-phylip-concat.charsets
sed -i 's/uce-/uce_/g' ./5-input-alignments/all-rev1/hirundinidae-all-rev1-95-nexus-concat/hirundinidae-all-rev1-95-nexus-concat.nexus
```

------------------------------------------------------------------------------------------
### 3. Full, revised 2 dataset (n = 103 taxa)

This includes the tissue samples and Pseudochelidon.

#### 2. Format concatenated alignment for full dataset

##### 1. Use phyluce to generate concatenated alignment phylip

```
cd /data3/hirundinidae_phylogeny/workflow_reanalysis/phyluce
phyluce_align_concatenate_alignments --alignments ./4-datasets/all-rev2/mafft-gblocks-nexus-internal-trimmed-all-rev2-samples-incomplete-clean-95p --output ./5-input-alignments/all-rev2/hirundinidae-all-rev2-95-phylip-concat --phylip --log-path .
```

##### 2. Use phyluce to generate concatenated alignment nexus

```
phyluce_align_concatenate_alignments --alignments ./4-datasets/all-rev2/mafft-gblocks-nexus-internal-trimmed-all-rev2-samples-incomplete-clean-95p --output ./5-input-alignments/all-rev2/hirundinidae-all-rev2-95-nexus-concat --nexus --log-path .
```

##### 3. Format UCE labels

```
sed -i 's/uce-/uce_/g' ./5-input-alignments/all-rev2/hirundinidae-all-rev2-95-phylip-concat/hirundinidae-all-rev2-95-phylip-concat.charsets
sed -i 's/uce-/uce_/g' ./5-input-alignments/all-rev2/hirundinidae-all-rev2-95-nexus-concat/hirundinidae-all-rev2-95-nexus-concat.nexus
```

## Part 7 - Maximum likelihood analysis in RAxML

With a formatted input dataset, we'll now estimate ML trees from the UCE data.

We have locus-specific alignments and a concatenated alignment as inputs.

### Set up environment

We'll make an analysis directory for RAxML and set up subdirectories for various analyses.
```
cd /data3/hirundinidae_phylogeny/workflow_renalysis/
mkdir raxml
```

------------------------------------------------------------------------------------------
### 1. Perform pilot concatenated analysis on all samples (n = 122 taxa)

This is based on 100 random UCE loci from the 95% UCE matrix.

The input data are in `/data3/hirundinidae_phylogeny/workflow_reanalysis/phyluce/5-input-alignments/all/hirundinidae-all-95-pilot-phylip-concat/hirundinidae-all-95-pilot-phylip-concat.phylip`

#### Set up environment

```
cd raxml
mkdir hirundinidae-all-95-concat-pilot
```

#### Retrieve the input concatenated data matrix

```
cd hirundinidae-all-95-concat-pilot
cp /data3/hirundinidae_phylogeny/workflow_reanalysis/phyluce/5-input-alignments/all/hirundinidae-all-95-pilot-phylip-concat/hirundinidae-all-95-pilot-phylip-concat.* .
mv hirundinidae-all-95-pilot-phylip-concat.phylip hirundinidae-all-95-pilot-concat.phylip
mv hirundinidae-all-95-pilot-phylip-concat.charsets hirundinidae-all-95-pilot-concat.charsets
```

#### 1. Perform initial analysis

We want to check to make sure RAxML reads in the data and will run.

```
raxml-ng-mpi --parse --msa ./hirundinidae-all-95-pilot-concat.phylip --model GTR+G
```

#### 2. Perform tree search

We'll do 40 ML tree searches (20 from a random starting tree and 20 from a parsimony starting tree). These steps follow Clare's notes in '9d-all-95' in her '9-RAxML-NG' folder.

```
raxml-ng-mpi --threads 36 --msa hirundinidae-all-95-pilot-concat.phylip.raxml.rba --model GTR+G --prefix ha95-pilot-best --seed 1487 --tree pars{20},rand{20}
```

#### 3. Bootstrapping

We'll run 100 bootstraps on the results of the tree search.

```
raxml-ng-mpi --threads 36 --bootstrap --msa hirundinidae-all-95-pilot-concat.phylip.raxml.rba --model GTR+G --prefix ha95-pilot-boot --seed 1913 --bs-trees 100
```

#### 4. Compute support values

First calculate bootstrap support on the best-scoring ML tree from RAxML.

```
raxml-ng-mpi --threads 1 --support --tree ha95-pilot-best.raxml.bestTree --bs-trees ha95-pilot-boot.raxml.bootstraps --prefix ha95-pilot-annotated
```

Then, summarize bootstraps as a 50 percent majority rule consensus tree using `sumtrees.py` from Jeet Sukumaran's [DendroPy code](https://github.com/jeetsukumaran/DendroPy/blob/main/).

Install DendroPy in python 3 environment:
```
conda activate pixy
pip install dendropy
```

Run SumTrees:
```
sumtrees.py -o ha95-pilot-consensus.sumtrees.support ha95-pilot-boot.raxml.bootstraps
```

------------------------------------------------------------------------------------------
### 2. Perform pilot concatenated analysis on the full, revised 1 dataset (n = 118 taxa)

This is based on all UCE loci from the 95% UCE matrix.

The input data are in `/data3/hirundinidae_phylogeny/workflow_reanalysis/phyluce/5-input-alignments/all-rev1/hirundinidae-all-rev1-95-pilot-phylip-concat/hirundinidae-all-rev1-95-pilot-phylip-concat.phylip`

#### Set up environment

```
cd raxml
mkdir hirundinidae-all-rev1-95-concat-pilot
```

#### Retrieve the input concatenated data matrix

```
cd hirundinidae-all-rev1-95-concat-pilot
cp /data3/hirundinidae_phylogeny/workflow_reanalysis/phyluce/5-input-alignments/all-rev1/hirundinidae-all-rev1-95-pilot-phylip-concat/hirundinidae-all-rev1-95-pilot-phylip-concat.* .
mv hirundinidae-all-rev1-95-pilot-phylip-concat.phylip hirundinidae-all-rev1-95-pilot-concat.phylip
mv hirundinidae-all-rev1-95-pilot-phylip-concat.charsets hirundinidae-all-rev1-95-pilot-concat.charsets
```

#### 1. Perform initial analysis

We want to check to make sure RAxML reads in the data and will run.

```
raxml-ng-mpi --parse --msa ./hirundinidae-all-rev1-95-pilot-concat.phylip --model GTR+G
```

#### 2. Perform tree search

We'll do 40 ML tree searches (20 from a random starting tree and 20 from a parsimony starting tree). These steps follow Clare's notes in '9d-all-95' in her '9-RAxML-NG' folder.

```
raxml-ng-mpi --threads 36 --msa hirundinidae-all-rev1-95-pilot-concat.phylip.raxml.rba --model GTR+G --prefix ha95-rev1-pilot-best --seed 1487 --tree pars{20},rand{20}
```

#### 3. Bootstrapping

We'll run 100 bootstraps on the results of the tree search.

```
raxml-ng-mpi --threads 36 --bootstrap --msa hirundinidae-all-rev1-95-pilot-concat.phylip.raxml.rba --model GTR+G --prefix ha95-rev1-pilot-boot --seed 2048 --bs-trees 100
```

#### 4. Compute support values

First calculate bootstrap support on the best-scoring ML tree from RAxML.

```
raxml-ng-mpi --threads 1 --support --tree ha95-rev1-pilot-best.raxml.bestTree --bs-trees ha95-rev1-pilot-boot.raxml.bootstraps --prefix ha95-rev1-pilot-annotated
```

Then, summarize bootstraps as a 50 percent majority rule consensus tree using `sumtrees.py` from Jeet Sukumaran's [DendroPy code](https://github.com/jeetsukumaran/DendroPy/blob/main/).

```
sumtrees.py -o ha95-rev1-pilot-consensus.sumtrees.support ha95-rev1-pilot-boot.raxml.bootstraps
```

------------------------------------------------------------------------------------------
### 3. Perform full concatenated analysis on the full, revised 1 dataset (n = 118 taxa)

Here we'll run RAxML on the concatenated 95% matrix for the all-rev1 dataset (n = 118 taxa; see above).

The input data are in `/data3/hirundinidae_phylogeny/workflow_reanalysis/phyluce/5-input-alignments/all-rev1/hirundinidae-all-rev1-95-phylip-concat/hirundinidae-all-rev1-95-phylip-concat.phylip`

#### Set up environment

```
cd raxml
mkdir hirundinidae-all-rev1-95-concat
```

#### Retrieve the input concatenated data matrix

```
cd hirundinidae-all-rev1-95-concat
cp /data3/hirundinidae_phylogeny/workflow_reanalysis/phyluce/5-input-alignments/all-rev1/hirundinidae-all-rev1-95-phylip-concat/hirundinidae-all-rev1-95-phylip-concat.* .
mv hirundinidae-all-rev1-95-phylip-concat.phylip hirundinidae-all-rev1-95-concat.phylip
mv hirundinidae-all-rev1-95-phylip-concat.charsets hirundinidae-all-rev1-95-concat.charsets
```

#### 1. Perform initial test run

We want to check to make sure RAxML reads in the data and will run.

```
raxml-ng-mpi --parse --msa ./hirundinidae-all-rev1-95-concat.phylip --model GTR+G
```

This completed and confirmed that RAxML could run on the input data. It produced the output files:
```
hirundinidae-all-rev1-95-concat.phylip.raxml.log
hirundinidae-all-rev1-95-concat.phylip.raxml.rba
```

#### 2. Perform tree search

We'll do 40 ML tree searches (20 from a random starting tree and 20 from a parsimony starting tree). These steps follow Clare's notes in '9d-all-95' in her '9-RAxML-NG' folder.

```
raxml-ng-mpi --threads 56 --msa hirundinidae-all-rev1-95-concat.phylip.raxml.rba --model GTR+G --prefix ha95-rev1-best --seed 1946 --tree pars{20},rand{20}
```

Note [08.09.2023]: Running this in mirrored directory on Nostromo to take advantage of extra threads.

#### 3. Bootstrapping

We'll run 100 bootstraps on the results of the tree search.

```
raxml-ng-mpi --threads 56 --bootstrap --msa hirundinidae-all-rev1-95-concat.phylip.raxml.rba --model GTR+G --prefix ha95-rev1-boot --seed 6985 --bs-trees 100
```

#### 4. Compute support values

First calculate bootstrap support on the best-scoring ML tree from RAxML.

```
raxml-ng-mpi --threads 1 --support --tree ha95-rev1-best.raxml.bestTree --bs-trees ha95-rev1-boot.raxml.bootstraps --prefix ha95-rev1-annotated
```

Then, summarize bootstraps as a 50 percent majority rule consensus tree using `sumtrees.py` from Jeet Sukumaran's [DendroPy code](https://github.com/jeetsukumaran/DendroPy/blob/main/).

```
conda activate pixy
sumtrees.py -o ha95-rev1-consensus.sumtrees.support ha95-rev1-boot.raxml.bootstraps
```

We'll incorporate the bootstrapped, best-scoring ML tree in the main results.

------------------------------------------------------------------------------------------
### 4. Perform full concatenated analysis on the full, revised 2 dataset (n = 103 taxa)

Here we'll run RAxML on the concatenated 95% matrix for the all-rev2 dataset (tissue samples and Pseudochelidon).

The input data are in `/data3/hirundinidae_phylogeny/workflow_reanalysis/phyluce/5-input-alignments/all-rev2/hirundinidae-all-rev2-95-phylip-concat/hirundinidae-all-rev2-95-phylip-concat.phylip`

#### Set up environment

```
cd raxml
mkdir hirundinidae-all-rev2-95-concat
```

#### Retrieve the input concatenated data matrix

```
cd hirundinidae-all-rev2-95-concat
cp /data3/hirundinidae_phylogeny/workflow_reanalysis/phyluce/5-input-alignments/all-rev2/hirundinidae-all-rev2-95-phylip-concat/hirundinidae-all-rev2-95-phylip-concat.* .
mv hirundinidae-all-rev2-95-phylip-concat.phylip hirundinidae-all-rev2-95-concat.phylip
mv hirundinidae-all-rev2-95-phylip-concat.charsets hirundinidae-all-rev2-95-concat.charsets
```

#### 1. Perform initial test run

Note: running this in mirrored directory on Nostromo.

We want to check to make sure RAxML reads in the data and will run.

```
raxml-ng-mpi --parse --msa ./hirundinidae-all-rev2-95-concat.phylip --model GTR+G
```

This completed and confirmed that RAxML could run on the input data. It produced the output files:
```
hirundinidae-all-rev2-95-concat.phylip.raxml.log
hirundinidae-all-rev2-95-concat.phylip.raxml.rba
```

#### 2. Perform tree search

We'll do 40 ML tree searches (20 from a random starting tree and 20 from a parsimony starting tree). These steps follow Clare's notes in '9d-all-95' in her '9-RAxML-NG' folder.

```
raxml-ng-mpi --threads 24 --msa hirundinidae-all-rev2-95-concat.phylip.raxml.rba --model GTR+G --prefix ha95-rev2-best --seed 1946 --tree pars{20},rand{20}
```

Note [08.09.2023]: Running this in mirrored directory on Nostromo to take advantage of extra threads.

#### 3. Bootstrapping

We'll run 100 bootstraps on the results of the tree search.

```
raxml-ng-mpi --threads 56 --bootstrap --msa hirundinidae-all-rev2-95-concat.phylip.raxml.rba --model GTR+G --prefix ha95-rev2-boot --seed 547689 --bs-trees 100
```

#### 4. Compute support values

First calculate bootstrap support on the best-scoring ML tree from RAxML.

```
raxml-ng-mpi --threads 1 --support --tree ha95-rev2-best.raxml.bestTree --bs-trees ha95-rev2-boot.raxml.bootstraps --prefix ha95-rev2-annotated
```

Then, summarize bootstraps as a 50 percent majority rule consensus tree using `sumtrees.py` from Jeet Sukumaran's [DendroPy code](https://github.com/jeetsukumaran/DendroPy/blob/main/).

```
conda activate pixy
sumtrees.py -o ha95-rev2-consensus.sumtrees.support ha95-rev2-boot.raxml.bootstraps
```

We'll incorporate the bootstrapped, best-scoring ML tree in the main results.

## Part 8 - Species tree inference

We'll run the coalescent-based quartets method SVDquartets on the concatenated alignment for the full, revised 1 dataset (n = 118 taxa).

With alignments for individual UCEs, we'll also estimate the species tree based on tissue samples using ASTRAL, first estimating gene trees using RAxML, then performing multi-species coalescent species tree inference using ASTRAL.

------------------------------------------------------------------------------------------
### Set up environment

```
cd /data3/hirundinidae_phylogeny/workflow_reanalysis/
mkdir species_tree
cd species_tree
mkdir svdq
mkdir astral
```

------------------------------------------------------------------------------------------
### Install SVDQuartets

SVDQuartets is implemented in PAUP*, so we'll install the whole package.

```
cd ~/tmp/
wget http://phylosolutions.com/paup-test/paup4a168_ubuntu64.gz
gunzip paup4a168_ubuntu64.gz
mv paup4a168_ubuntu64 paup
chmod +x paup
cp paup /usr/local/bin
```

------------------------------------------------------------------------------------------
### Install raxmlHPC-AVX

#### Create new conda environment and install RAxML (separate from raxml-ng in base environment)

```
conda create --name raxml
conda activate raxml
conda install -c bioconda raxml
```

This installs `raxmlHPC-AVX2`, which has the same input commands as the previous version used by Clare.

Note, with `raxmlHPC-PTHREADS-AVX2`, one can also specify the `-T` flag with an integer number of threads for parallel processing.

------------------------------------------------------------------------------------------
### Install Astral

Astral is available from https://github.com/smirarab/ASTRAL.

#### Set up environment and install Astral

```
cd astral
git clone https://github.com/smirarab/ASTRAL.git
cd ASTRAL
./make.sh
```

#### Test installation

```
cd Astral
java -jar astral.5.7.8.jar -i ./main/test_data/song_primates.424.gene.tre
```

It works.

------------------------------------------------------------------------------------------
### 1. Pilot species tree analysis (SVDquartets; all-rev1 dataset, n = 118 taxa, 100 loci)

Brant recommended SVDquartets as our main coalescent-based inference to compare to the concatenated analysis.

The ASTRAL analysis below can be used on tissue samples only, and mainly to evaluate concordance in higher-level relationships.

#### Set up environment

```
cd /data3/hirundinidae_phylogeny/workflow_reanalysis/species_tree/
mkdir svdq
cd svdq
mkdir ha95-rev1-pilot
```

#### Retrieve the input concatenated data matrix

```
cd ha95-rev1-pilot
cp /data3/hirundinidae_phylogeny/workflow_reanalysis/phyluce/5-input-alignments/all-rev1/hirundinidae-all-rev1-95-pilot-nexus-concat/hirundinidae-all-rev1-95-pilot-nexus-concat.nexus .
mv hirundinidae-all-rev1-95-pilot-nexus-concat.nexus hirundinidae-all-rev1-95-pilot-concat.nexus
```

#### Run preliminary analysis, specifying outgroup clade

We'll run the analysis with 100 bootstrap replicates.
```
paup
exe hirundinidae-all-rev1-95-pilot-concat.nexus;
outgroup 25;
outgroup 42;
set outroot=mono;
svdq
svdq showScores=no seed=1234568 bootstrap nreps=100 treeFile=ha95-rev1-pilot.bootstrap.tre;
```

Save the tree to file:
```
savetrees file=ha95-rev1-pilot.svdq.tre savebootp=nodelabels;
```

------------------------------------------------------------------------------------------
### 2. Species tree analysis (SVDquartets; all-rev1 dataset, n = 118 taxa, all loci)

This is the main coalescent-based species tree analysis to report on.

#### Set up environment

```
cd /data3/hirundinidae_phylogeny/workflow_reanalysis/species_tree/svdq
mkdir ha95-rev1
```

#### Retrieve the input concatenated data matrix

```
cd ha95-rev1
cp /data3/hirundinidae_phylogeny/workflow_reanalysis/phyluce/5-input-alignments/all-rev1/hirundinidae-all-rev1-95-nexus-concat/hirundinidae-all-rev1-95-nexus-concat.nexus .
mv hirundinidae-all-rev1-95-nexus-concat.nexus hirundinidae-all-rev1-95-concat.nexus
```

#### Run preliminary analysis, specifying outgroup clade

We'll run the analysis with 100 bootstrap replicates.
```
paup
exe hirundinidae-all-rev1-95-concat.nexus;
outgroup 25;
outgroup 42;
set outroot=mono;
svdq showScores=no seed=65478211 bootstrap nreps=100 treeFile=ha95-rev1.bootstrap.tre;
```

Save the tree to file:
```
savetrees file=ha95-rev1.svdq.tre savebootp=nodelabels;
```

------------------------------------------------------------------------------------------
### 3. Species tree analysis (SVDquartets; all-rev2 dataset, tissue samples + Pseudochelidon)

This is for comparison to the all-rev1 dataset results.

#### Set up environment

```
cd /data3/hirundinidae_phylogeny/workflow_reanalysis/species_tree/svdq
mkdir ha95-rev2
```

#### Retrieve the input concatenated data matrix

```
cd ha95-rev2
cp /data3/hirundinidae_phylogeny/workflow_reanalysis/phyluce/5-input-alignments/all-rev2/hirundinidae-all-rev2-95-nexus-concat/hirundinidae-all-rev2-95-nexus-concat.nexus .
mv hirundinidae-all-rev2-95-nexus-concat.nexus hirundinidae-all-rev2-95-concat.nexus
```

#### Run preliminary analysis, specifying outgroup clade

We'll run the analysis with 100 bootstrap replicates.
```
paup
exe hirundinidae-all-rev2-95-concat.nexus;
outgroup 23;
outgroup 37;
set outroot=mono;
svdq showScores=no seed=65478211 bootstrap nreps=100 treeFile=ha95-rev2.bootstrap.tre;
```

Save the tree to file:
```
savetrees file=ha95-rev2.svdq.tre savebootp=nodelabels;
```

------------------------------------------------------------------------------------------
### 4. Species tree analysis (ASTRAL; all-rev3 dataset; tissue samples only)

This includes the tissue dataset (n = 101 taxa). Preliminary analysis revealed that Pseudochelidon gets placed within the clade with Hirundo, Delichon, etc., which can't be real and is probably a data artifact.

Paring down to just the tissue samples gives us the ability to compare the topologies for Hirundininae across approaches (i.e., no Pseudochelidon included).

#### Set up environment

```
cd /data3/hirundinidae_phylogeny/workflow_reanalysis/species_tree/astral
mkdir ha95-rev3
cd ha95-rev3
mkdir raxml
cd raxml
mkdir phylip
mkdir gene_trees
```

#### 1. Gene tree estimation

We'll estimated gene trees using RAxML, which will be used as input for the coalescent-based method in Astral.

##### 1. Convert input alignments

We first need to convert nexus alignments from phyluce to relaxed phylip input for RAxML.

We'll run the conversion in the `phyluce` directory, then copy over to our working directory for gene tree estimation.

```
cd /data3/hirundinidae_phylogeny/workflow_reanalysis/phyluce/
conda activate phyluce
phyluce_align_convert_one_align_to_another --cores 8 --alignments ./4-datasets/all-rev3/mafft-gblocks-nexus-internal-trimmed-all-rev3-samples-incomplete-clean-95p --output ./5-input-alignments/all-rev3/hirundinidae-all-rev3-95-phylip --input-format nexus --output-format phylip-relaxed --log-path .
cd /data3/hirundinidae_phylogeny/workflow_reanalysis/species_tree/astral/ha95-rev3/raxml/phylip
cp /data3/hirundinidae_phylogeny/workflow_reanalysis/phyluce/5-input-alignments/all-rev3/hirundinidae-all-rev3-95-phylip/*.phylip-relaxed .
cd ..
```

##### 2. Count bootstrap replicates needed for multilocus bootstrapping

We will generate consensus trees, and rather than running bootstraps on each tree, we'll perform multilocus boostrapping.

We need to know how many bootstrap replicates we'll need per locus, and can figure this out using `phyluce_genetrees_generate_multilocus_bootstrap_count`.

```
phyluce_genetrees_generate_multilocus_bootstrap_count --alignments ./phylip --bootstrap_replicates ha95-rev3.bootstrap.replicates --bootstrap_counts ha95-rev3.bootstrap.counts --directory ./phylip --bootreps 500 --log-path .  
```

This outputs two output files, one with the number of replicates needed for each locus and the other with bootstrap replicate sampling:
```
ha95-rev3.bootstrap.counts
ha95-rev3.bootstrap.replicates
```

##### 3. Perform test analysis

```
raxmlHPC-PTHREADS-AVX2 -T 48 -m GTRGAMMA -N 523 -p 125792 -b 399802 -n bootrep -k -s ./phylip/uce-7236.phylip-relaxed
```

This outputs the files:
```
RAxML_bootstrap.bootrep
RAxML_info.bootrep
```

These files are written to the current directory, so in the script below, we'll temporarily rename these and move to locus-specific output directories.

##### 4. Perform gene tree bootstrapping analysis

###### Set up environment

```
cd gene_trees
mkdir bootstrap
cd ..
```

###### Wrote `runRAxML_bootstrap_parallel.sh` to perform analysis on each UCE locus

This uses GNU parallel to run 48 simultaneous analyses.

The results are written to `./gene_trees/bootstrap/`.

runRAxML_bootstrap_parallel.sh:

```
file=$1
reps=$2
rand1=`shuf -i 1-1000000 -n 1`
rand2=`shuf -i 1-1000000 -n 1`
uce=`echo $file | cut -d'/' -f3 | cut -d'.' -f1`
mkdir ./gene_trees/bootstrap/$uce
raxmlHPC-AVX2 -m GTRGAMMA -N $reps -p $rand1 -b $rand2 -n $uce.bootrep -k -s $file
mv ./RAxML_bootstrap.$uce.bootrep ./gene_trees/bootstrap/$uce/RAxML_bootstrap.bootrep
mv ./RAxML_info.$uce.bootrep ./gene_trees/bootstrap/$uce/RAxML_info.bootrep
```

###### Run the script with parallel

This reads in the bootstrap count output from phyluce as input.

```
conda activate raxml
parallel --progress --joblog logfile.raxml -j 48 --workdir . --colsep ' ' ./runRAxML_bootstrap_parallel.sh :::: ha95-rev3.bootstrap.counts
```

###### Sort bootstrap replicates

We'll use `phyluce_genetrees_sort_multilocus_bootstraps` to sort bootstrapped gene trees.

```
cd gene_trees
mkdir bootstrap_sorted
cd ..
phyluce_genetrees_sort_multilocus_bootstraps --input ./gene_trees/bootstrap --bootstrap_replicates ./ha95-rev3.bootstrap.replicates --output ./gene_trees/bootstrap_sorted
```

##### 5. Perform coalescent-based species tree analysis in Astral

###### Set up environment

Get into the right place and retrieve bootstrap trees:
```
cd /data3/hirundinidae_phylogeny/workflow_reanalysis/species_tree/astral/ha95-rev3
mkdir astral_trees
mkdir astral_log
```

Make list file with path to multilocus bootstraps:
```
cd ./raxml/gene_trees
list=`ls -lh bootstrap_sorted/boot* | cut -d':' -f2 | cut -d' ' -f2`; for i in $list; do echo /data3/hirundinidae_phylogeny/workflow_reanalysis/species_tree/astral/ha95-rev3/raxml/gene_trees/$i >> ../../bootstrap.list; done
```

Note: change the path to `/data3/` instead of `/media/mother/extradrive1` if working on Terminator.

###### Run test analysis

```
java -Xmx3900M -jar ../ASTRAL/astral.5.7.8.jar -i /media/mother/extradrive1/hirundinidae_phylogeny/workflow_reanalysis/species_tree/astral/ha95-rev3/raxml/gene_trees/bootstrap_sorted/boot000 -o test.boot000.astral.tre
```

###### Write run script for parallel analysis

runAstral_parallel.sh:

```
file=$1
name=`echo $file | rev | cut -d"/" -f1 | rev`
java -Xmx3900M -jar ../ASTRAL/astral.5.7.8.jar -i $file -o ./astral_trees/$name.astral.tre 2>./astral_log/$name.astral.log
```

###### Run Astral on pilot dataset

```
parallel --progress --joblog logfile.astral -j 36 --workdir . ./runAstral_parallel.sh :::: bootstrap.list
```

###### Concatenate Astral bootrep species trees

```
cat ./astral_trees/*.tre > astral.ha95-rev3.bootrep.tre
```

###### Summarize trees using 50% majority rule consensus tree

```
conda activate pixy
sumtrees.py -o astral.ha95-rev3.50con.tre --summary-target=consensus --min-clade-freq=0.5 astral.ha95-rev3.bootrep.tre
```

## Part 9 - Divergence dating using Bayesian inference

We'll use beast2 to perform full Bayesian tree search and divergence dating analyses on subsets of loci.

Useful resources:

https://beast2-dev.github.io/beast-docs/beast2/DivergenceDating/DivergenceDatingTutorial.html
https://github.com/jasonleebrown/UCE_phyluce_pipeline/blob/master/README.md#bayesian-analysis-with-beast
https://beast2-dev.github.io/hmc/hmc//Standard/Clock_Model/index.html
https://beast.community/errors.html

------------------------------------------------------------------------------------------
### Set up environment

```
cd /data3/hirundinidae_phylogeny/workflow_renalysis/
mkdir beast
cd beast
mkdir ha95-rev1
```

------------------------------------------------------------------------------------------
### Install beast2

```
cd ~/tmp/
wget https://github.com/CompEvol/beast2/releases/download/v2.7.5/BEAST.v2.7.5.Linux.x86.tgz
tar -xf BEAST.v2.7.5.Linux.x86.tgz
mkdir /data3/hirundinidae_phylogeny/workflow_reanalysis/beast
cp -r ./beast /data3/hirundinidae_phylogeny/workflow_reanalysis/beast/
```

Also run included `packagemanager` to install ORC package:
```
./beast/bin/packagemanager -add ORC
```

------------------------------------------------------------------------------------------
### Install Beagle library

```
cd ~/tmp/
git clone --depth=1 https://github.com/beagle-dev/beagle-lib.git
cd beagle-lib
mkdir build
cd build
cmake -DCMAKE_INSTALL_PREFIX:PATH=$HOME ..
make install
```

This installed the Beagle libraries to `~/lib/`.

------------------------------------------------------------------------------------------
### 1. Select random sets of UCE loci for analysis

We'll select 5 sets of 50 random loci.

```
cd /data3/hirundinidae_phylogeny/workflow_renalysis/phyluce
mkdir ./4-datasets/all-rev1/mafft-gblocks-nexus-internal-trimmed-all-rev1-samples-incomplete-clean-95p-subset
cd ./4-datasets/all-rev1/mafft-gblocks-nexus-internal-trimmed-all-rev1-samples-incomplete-clean-95p-subset
mkdir subset1 subset2 subset3 subset4 subset5
cd ../mafft-gblocks-nexus-internal-trimmed-all-rev1-samples-incomplete-clean-95p
list=`ls -lh *.nexus | cut -d':' -f2 | cut -d' ' -f2 | shuf -n 50`; for i in $list; do cp $i ../mafft-gblocks-nexus-internal-trimmed-all-rev1-samples-incomplete-clean-95p-subset/subset1; done
list=`ls -lh *.nexus | cut -d':' -f2 | cut -d' ' -f2 | shuf -n 50`; for i in $list; do cp $i ../mafft-gblocks-nexus-internal-trimmed-all-rev1-samples-incomplete-clean-95p-subset/subset2; done
list=`ls -lh *.nexus | cut -d':' -f2 | cut -d' ' -f2 | shuf -n 50`; for i in $list; do cp $i ../mafft-gblocks-nexus-internal-trimmed-all-rev1-samples-incomplete-clean-95p-subset/subset3; done
list=`ls -lh *.nexus | cut -d':' -f2 | cut -d' ' -f2 | shuf -n 50`; for i in $list; do cp $i ../mafft-gblocks-nexus-internal-trimmed-all-rev1-samples-incomplete-clean-95p-subset/subset4; done
list=`ls -lh *.nexus | cut -d':' -f2 | cut -d' ' -f2 | shuf -n 50`; for i in $list; do cp $i ../mafft-gblocks-nexus-internal-trimmed-all-rev1-samples-incomplete-clean-95p-subset/subset5; done
```
 
------------------------------------------------------------------------------------------
### 2. Format input nexus files for locus subsets

#### 1. Generate concatenated alignment nexus files for each subset

```
phyluce_align_concatenate_alignments --alignments ./4-datasets/all-rev1/mafft-gblocks-nexus-internal-trimmed-all-rev1-samples-incomplete-clean-95p-subset/subset1/ --output ./5-input-alignments/all-rev1/hirundinidae-all-rev1-95-nexus-subset1 --nexus --log-path .
phyluce_align_concatenate_alignments --alignments ./4-datasets/all-rev1/mafft-gblocks-nexus-internal-trimmed-all-rev1-samples-incomplete-clean-95p-subset/subset2/ --output ./5-input-alignments/all-rev1/hirundinidae-all-rev1-95-nexus-subset2 --nexus --log-path .
phyluce_align_concatenate_alignments --alignments ./4-datasets/all-rev1/mafft-gblocks-nexus-internal-trimmed-all-rev1-samples-incomplete-clean-95p-subset/subset3/ --output ./5-input-alignments/all-rev1/hirundinidae-all-rev1-95-nexus-subset3 --nexus --log-path .
phyluce_align_concatenate_alignments --alignments ./4-datasets/all-rev1/mafft-gblocks-nexus-internal-trimmed-all-rev1-samples-incomplete-clean-95p-subset/subset4/ --output ./5-input-alignments/all-rev1/hirundinidae-all-rev1-95-nexus-subset4 --nexus --log-path .
phyluce_align_concatenate_alignments --alignments ./4-datasets/all-rev1/mafft-gblocks-nexus-internal-trimmed-all-rev1-samples-incomplete-clean-95p-subset/subset5/ --output ./5-input-alignments/all-rev1/hirundinidae-all-rev1-95-nexus-subset5 --nexus --log-path .
```

#### 2. Format UCE labels

```
sed -i 's/uce-/uce_/g' ./5-input-alignments/all-rev1/hirundinidae-all-rev1-95-nexus-subset1/hirundinidae-all-rev1-95-nexus-subset1.nexus
sed -i 's/uce-/uce_/g' ./5-input-alignments/all-rev1/hirundinidae-all-rev1-95-nexus-subset2/hirundinidae-all-rev1-95-nexus-subset2.nexus
sed -i 's/uce-/uce_/g' ./5-input-alignments/all-rev1/hirundinidae-all-rev1-95-nexus-subset3/hirundinidae-all-rev1-95-nexus-subset3.nexus
sed -i 's/uce-/uce_/g' ./5-input-alignments/all-rev1/hirundinidae-all-rev1-95-nexus-subset4/hirundinidae-all-rev1-95-nexus-subset4.nexus
sed -i 's/uce-/uce_/g' ./5-input-alignments/all-rev1/hirundinidae-all-rev1-95-nexus-subset5/hirundinidae-all-rev1-95-nexus-subset5.nexus
```

#### 3. Remove single quotes and '.nexus' from partition character sets

BEAUTi barfs on these and does not read in the nexus file if they are included.

```
sed -i "s/'//g" ./5-input-alignments/all-rev1/hirundinidae-all-rev1-95-nexus-subset1/hirundinidae-all-rev1-95-nexus-subset1.nexus
sed -i "s/'//g" ./5-input-alignments/all-rev1/hirundinidae-all-rev1-95-nexus-subset2/hirundinidae-all-rev1-95-nexus-subset2.nexus
sed -i "s/'//g" ./5-input-alignments/all-rev1/hirundinidae-all-rev1-95-nexus-subset3/hirundinidae-all-rev1-95-nexus-subset3.nexus
sed -i "s/'//g" ./5-input-alignments/all-rev1/hirundinidae-all-rev1-95-nexus-subset4/hirundinidae-all-rev1-95-nexus-subset4.nexus
sed -i "s/'//g" ./5-input-alignments/all-rev1/hirundinidae-all-rev1-95-nexus-subset5/hirundinidae-all-rev1-95-nexus-subset5.nexus
```

```
sed -i 's/.nexus//g' ./5-input-alignments/all-rev1/hirundinidae-all-rev1-95-nexus-subset1/hirundinidae-all-rev1-95-nexus-subset1.nexus
sed -i 's/.nexus//g' ./5-input-alignments/all-rev1/hirundinidae-all-rev1-95-nexus-subset2/hirundinidae-all-rev1-95-nexus-subset2.nexus
sed -i 's/.nexus//g' ./5-input-alignments/all-rev1/hirundinidae-all-rev1-95-nexus-subset3/hirundinidae-all-rev1-95-nexus-subset3.nexus
sed -i 's/.nexus//g' ./5-input-alignments/all-rev1/hirundinidae-all-rev1-95-nexus-subset4/hirundinidae-all-rev1-95-nexus-subset4.nexus
sed -i 's/.nexus//g' ./5-input-alignments/all-rev1/hirundinidae-all-rev1-95-nexus-subset5/hirundinidae-all-rev1-95-nexus-subset5.nexus
```

------------------------------------------------------------------------------------------
### 3. Prepare XML files in BEAUTi

#### Edit nexus alignments to contain single partitions

Do this manually based on the start of the first UCE partition and end of the last UCE partition.

Also remove the 'charpartition combined' line. The charsets block at the end should look something like:
```
begin sets;
charset uce = 1-48936;
end;
```

#### Set parameters in BEAUTi

Site model:
* Gamma site model
	* Substitution rate = 1
	* Gamma category count = 4
	* Shape = 1 (estimate)
	* Proportion invariant = 0
	* Substitution model = HKY
		* Kappa = 2 (estimate)
		* Frequencies = empirical
		
Clock model:
* Optimized relaxed clock
	* Mean clock rate = 1
	
Priors
* Defaults
* Add MRCA prior ('HDA' for Hirundo, Delichon, and Riparia)
	* Set as monophyletic
	* Include all taxa except for outgroups, Pseudochelidon, and Psalidoprocne
	* Distribution = normal
		* Mean = 9
		* Sigma = 1.5
		* Offset = 0
* Add MRCA prior ('ANC' for all taxa)
	* Set as monophyletic
	* Include all outgroup and ingroup taxa
	* Distribution = normal
		* Mean = 22
		* Sigma = 2.5
		* Offset = 0
		
MCMC
* Chain length = 10000000
* Store every = 1000
* Trace log every = 1000

#### Get XML files into correct analysis directory on Terminator

`/data3/hirundinidae_phylogeny/workflow_reanalysis/beast/ha95-rev1/`

------------------------------------------------------------------------------------------
### 4. Run BEAST analysis

We'll need to export the location of the Beagle library files, installed above.

```
cd /data3/hirundinidae_phylogeny/workflow_reanalysis/beast/ha95-rev1/
export LD_LIBRARY_PATH=~/lib:$LD_LIBRARY_PATH
nohup ../beast/bin/beast -threads 8 hirundinidae-all-rev1-95-nexus-subset1.xml > runlog.ha95-rev1-subset1.log &
nohup ../beast/bin/beast -threads 8 hirundinidae-all-rev1-95-nexus-subset2.xml > runlog.ha95-rev1-subset2.log &
nohup ../beast/bin/beast -threads 8 hirundinidae-all-rev1-95-nexus-subset3.xml > runlog.ha95-rev1-subset3.log &
nohup ../beast/bin/beast -threads 8 hirundinidae-all-rev1-95-nexus-subset4.xml > runlog.ha95-rev1-subset4.log &
nohup ../beast/bin/beast -threads 8 hirundinidae-all-rev1-95-nexus-subset5.xml > runlog.ha95-rev1-subset5.log &
```

## Part 10 - Historical Biogeography

We'll use BioGeoBEARS to test and compare biogeographic models to estimate the historical biogeography across the Hirundinidae tree.

The Bayesian time-calibrated phylogeny will be used as the topology for these analyses (pruned to n = 88 representative taxa).

### Set up environment

In Google Drive project folder
```
cd workflow_reanalysis/analysis
mkdir biogeography
cd biogeography
mkdir data
```

### 1. Format input tree topology

The time-calibrated tree was estimated using BEAST, with branch lengths in time (millions of years).

Format taxon short names for pruned n = 88 dataset in `./hirundinidae-all-rev1-95-prune.short-taxon-names.txt`.

Run `../../R/BioGeoBEARS_format_input_tree.R` to format input BEAST nexus tree to tree output with shortened taxon names.

### 2. Format input geographic data

The geographic data collected by Clare are in `swallow_data_species.xlsx`; copy it into `./data/`.

Format `./data/swallow_geographic_data_8realms.xlsx` to configure input data.
* This table contains rows for the 88 taxa, with presence/absence (1/0) for the 8 geographic realms used in the study.
* The columns need to be ordered as follows:
	* A = afrotropical
	* B = paleartic
	* C = oriental
	* D = australian
	* E = oceanian
	* F = nearctic
	* G = panamanian
	* H = neotropical

Then, format `/data/swallow_geographic_data_8realms.txt`; a simple tab-delimited text version of the table.

From this, format an input geographic data phylip file for BioGeoBEARs, `./data/swallow_geographic_data_8realms.phy`.

The first 7 lines:
```
88	8	(A B C D E F G H)
Alopochelidon_fucata	00000001
Atticora_fasciata	00000001
Atticora_pileata	00000010
Atticora_tibialis	00000011
Cecropis_abyssinica	10000000
Cecropis_cucullata	10000000
```

### 3. Format manual dispersal multiplier file

This file (`./data/swallow_manual_dispersal_multipliers.txt`) specifies how dispersal from/to different biogeographic realms can occur in the tested models.

The organization of this file is why the realm columns in the above files need to be in a particular order.

File contents:
```
A	B	C	D	E	F	G	H
1	1	1	0	0	0	0	0
1	1	1	1	1	1	0	0
1	1	1	1	1	0	0	0
0	1	1	1	1	0	0	0
0	1	1	1	1	0	0	0
0	1	0	0	0	1	1	1
0	0	0	0	0	1	1	1
0	0	0	0	0	1	1	1

END
```

### 4. BioGeoBEARS analysis

Run `../../R/BioGeoBEARS_constrained_dispersal.R` to import tree/geographic data, specify/run models, then perform model comparisons.

Running all 6 models will take ~1 hour.

The runs and data are saved to .Rdata files in the working directory, and the model comparison results are saved to various files:
```
Hirundinidae_constrained_dispersal_restable_AIC_rellike_formatted.txt
Hirundinidae_constrained_dispersal_restable_AIC_rellike.txt
Hirundinidae_constrained_dispersal_restable_AICc_rellike_formatted.txt
Hirundinidae_constrained_dispersal_restable_AICc_rellike.txt
Hirundinidae_constrained_dispersal_restable.txt
Hirundinidae_constrained_dispersal_teststable.txt
```

### 5. Plotting BioGeoBEARs results

Run `../../R/BioGeoBEARS_results_figure.R` to plot the results.

The script will point you to run `../../R/BioGeoBEARS_get_MLstates.R` to gather information for plotting ML states onto the tree nodes.

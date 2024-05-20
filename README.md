# Hirundinidae Phylogeny

![Consensus tree](/tree.png "cover image")

This repository contains details on the data processing and analysis steps used to estimate the swallow family phylogeny, divergence times, historical biogeography, and ancestral state reconstructions for sociality and nesting type. This workflow is a companion to the methods described in Schield & Brown et al. (_in Press_).

Note: this repository assumes a specific file organization. For reproducibility, your environment may vary and will need to be adjusted accordingly.

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

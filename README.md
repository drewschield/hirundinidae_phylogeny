# Hirundinidae Phylogeny

![Consensus tree](/tree.png "cover image")

This repository contains details on the data processing and analysis steps used to estimate the swallow family phylogeny, divergence times, historical biogeography, and ancestral state reconstructions for sociality and nesting type. This workflow is a companion to the methods described in Schield & Brown et al. (_in Press_).

Note: this repository assumes a specific file organization. For reproducibility, your environment may vary and will need to be adjusted accordingly.

# Workflow for analysis of UCE phylogeny for Hirundinidae - Part 0 - General information and installation of phyluce

The dataset for this study consists of UCE loci sampled from a majority of swallow species.

The samples include preserved tissue and historical skins, the latter which yielded lower quality genomic data, resulting in spurious results and very long branch lengths for some samples from skins (and several that were dropped from analysis completely).

------------------------------------------------------------------------------------------
### Set up environment

```
cd /data3/hirundinidae_phylogeny
mkdir workflow_reanalysis
```

------------------------------------------------------------------------------------------
## Install phyluce

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
## Install CIAlign

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

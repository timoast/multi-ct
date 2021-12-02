# multiCUT&Tag

Re-analysis of multiCUT&Tag data from https://doi.org/10.1016/j.molcel.2021.09.019

Install dependencies

```
mamba env create -f env.yaml
conda activate mct
```

This pipeline involves downloading data from AWS, which includes download costs. Make sure that you have an AWS account configured.

```
snakemake -j 24
```

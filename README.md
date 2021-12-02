# multiCUT&Tag

Re-analysis of single-cell multiCUT&Tag data from Gopalan et al., _Molecular Cell_ (2021).

Paper: https://doi.org/10.1016/j.molcel.2021.09.019  
GEO: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE171554  
SRA: https://www.ncbi.nlm.nih.gov/sra?term=SRP313652  
Code: https://github.com/snehagopalan710/Bulk-Multi-CUT-Tag

Install dependencies

```
mamba env create -f env.yaml
conda activate mct
```

This pipeline involves downloading data from AWS, which includes download costs. Make sure that you have an AWS account configured.

```
snakemake -j 24
```

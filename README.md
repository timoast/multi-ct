# multiCUT&Tag

Re-analysis of single-cell multiCUT&Tag data from Gopalan et al., _Molecular Cell_ (2021).

Paper: https://doi.org/10.1016/j.molcel.2021.09.019  
GEO: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM5227096  
SRA: https://trace.ncbi.nlm.nih.gov/Traces/sra/?run=SRR14150925  
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

## Steps

1. Download raw data
2. Attach cell barcodes
3. Demultiplex based on Tn5 barcodes and trim reads
4. Map reads to mm10 genome
5. Create fragment file

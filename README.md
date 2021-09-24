# radseq_snakemake
Snakemake scaffold for RADSEQ

## Getting started

```
conda env create --name radseq --file environment.yml
conda activate radseq

snakemake -j 16 --use-conda --rerun-incomplete --latency-wait 15 --resources mem_mb=200000 --cluster "sbatch -t 10080 -J radseq -p bmm -n 1 -N 1 -c {threads} --mem={resources.mem_mb}" -k -n
```

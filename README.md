# skip - sola klamydia illumina pipeline

### Install

1. `git clone` this repo:

```
git clone https://github.com/ammaraziz/skip
```

2. Install `miniconda` or preferrably `mamba`
	- https://mamba.readthedocs.io/en/latest/installation.html
	- https://docs.conda.io/projects/miniconda/en/latest/

3. Install `snakemake`:
```
mamba install -c bioconda snakemake 
```

4. Download the human genome, rename to `resources/GRCh38.fasta`

5. Done - The pipeline will handle the dependencies internally.

### Usage

1. Create a folder `skip/input/`
2. Add your fastq.gz files in `skip/input/`. Ensure they're named as follows: `{sample_name}_{direction}.fastq.gz`. eg 
`SRR12345_R1.fastq.gz` and `SRR12345_R2.fastq.gz`.
3. In `skip/` folder run the pipeline:
```
snakemake -j 8 --use-conda -k
```

- `-j 8` specifies the number of threads to use in total.
- `--use-conda` tells snakemake to install the dependencies.
- `-k` tells snakemake to keep going if a sample fails.

### Cite

Pipeline is created by Shola Olagoke with assistance from Ammar Aziz.

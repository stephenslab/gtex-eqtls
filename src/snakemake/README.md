# Multi-tissue EQTL Workflow
## Create Environment with Required Software
Commands below creates the `Python`/`R` environment required for this project

```bash
  conda create -n mt-eqtl -c bioconda --file python-packages.txt
  source activate mt-eqtl # to deactivate use: source deactivate
  conda install -c r r=3.2.2
  Rscript R-packages.R
```

On Midway cluster several modules should be loaded

```bash
  module load libtool/2.4
  module load zlib/1.2
  module load gsl/1.16
```

### Dependencies
```
  python=3.5.0
  snakemake=3.4.2
  scipy=0.16.0
  pandas=0.17.0
  h5py=2.5.0
  pytables=3.2.2
```

```r
  repo <- "http://cran.us.r-project.org"
  install.packages('devtools', repos = repo)
  require(devtools)
  install_version("ggplot2", version = "1.0.1", repos = repo)
  install_version("RSQLite", version = "1.0.0", repos = repo)
  source("http://bioconductor.org/biocLite.R")
  biocLite("rhdf5") # too bad I cannot install specified version for rhdf5
```

## The Pipeline
Implementations are in `rules` and `workflow` directories listed below. The usage are found in project [documentation](../../doc/notes)

```
  ../../src/snakemake
  |-- config.bashrc
  |-- python-packages.txt
  |-- README.md
  |-- R-packages.R
  |-- rules
  |   |-- eqtlbma-bf.rules
  |   |-- file-processing.rules
  |   |-- matrix-eqtl.rules
  |   `-- sumstat-to-h5.rules
  `-- workflows
      |-- cluster.yaml
      |-- config.yaml -> midway-V6-lite.yaml
      |-- eqtlbma
      |   `-- Snakefile
      |-- Makefile
      |-- midway-V6-full.yaml
      |-- midway-V6-lite.yaml
      `-- preprocessing
          `-- Snakefile
  
  4 directories, 15 files
  
```

Available rules are

```
  snakemake sumstat_to_h5
  snakemake prepare_merge_batch
  snakemake merge_h5
  snakemake count_genes_from_data
  snakemake count_genes_from_h5
  snakemake sample_max_null
  snakemake prepare_snp_lookup_db
  snakemake create_snp_lookup_db
  snakemake prepare_matrix_eqtl_input
  snakemake find_genes_with_sumstats
  snakemake eqtlbma_toy
  snakemake extract_data
  snakemake prepare_coords
  snakemake prepare_input_lists
  snakemake eqtlbma_batch_poisson
  
```

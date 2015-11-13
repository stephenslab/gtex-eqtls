# Multi-tissue EQTL Workflow
## Create Environment with Required Software
Commands below creates the `Python`/`R` environment required for this project

```bash
  conda create -n mt-eqtl -c bioconda --file python-packages.txt
  source activate mt-eqtl # to deactivate use: source deactivate
  conda install -c r r=3.2.2
  Rscript R-packages.R
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

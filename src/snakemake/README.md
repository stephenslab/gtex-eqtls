# Multi-tissue EQTL Workflow
## Create Environment with Required Software
Commands below creates the `Python`/`R` environment required for this project

```bash
  conda create -n mt-eqtl -c bioconda --file python-packages.txt
  source activate mt-eqtl # to deactivate use: source deactivate
  conda install -c r r=3.2.2
  Rscript R-packages.R
```

# Welcome to GTEx V8 multivariate analysis!

## Mashing the GTEx V8 data

For those interested in reproducing the MASH analysis on GTEx V8 data, 
[this notebook][mashr_wf] contains (should contain) sufficient information to get the analysis done.
Computational steps in the notebook are also outlined below. 
As you can see from the embedded links, some steps are workflows written and executed from other notebooks.

1. [Convert `fastqtl` summary statistics to `mashr` input][fastqtl_wf]
2. [Fit MASH model][mashr_wf]
3. [Prepare a subset of gene-SNP pairs to follow up with `mashr`][indep_qtl_wf]
4. [Compute posterior under MASH model][mashr_wf]

## A note to repo contributors

Type the following command to build this website from `Rmd` / `ipynb` files:

```bash
./release
```

and checkout the results here at:

https://gaow.github.io/mnm-gtex-v8

[mashr_wf]: https://gaow.github.io/mnm-gtex-v8/analysis/mashr_flashr_workflow.html
[fastqtl_wf]: https://gaow.github.io/mnm-gtex-v8/analysis/fastqtl_to_mash.html
[indep_qtl_wf]: https://gaow.github.io/mnm-gtex-v8/analysis/Independent_eQTL_Results.html

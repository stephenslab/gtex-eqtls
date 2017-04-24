# General Information
## Release & Analysis V4
**_Obsolete_**

*	 GTEx raw and normalized expression and PEER factors data:
	*	 ftp://ftp.broadinstitute.org/GTEx_Analysis_2014-06-13/expressionCovariateSandbox/
	*	 There is a `README` file in this folder that contains description of data
*	 GTEx summary statistics formatted to eqtlbma input, by Sarah:
	*	 [Midway] /project/mstephens/gtex/preprocessing/june2014/sumstats
*	 GTEx summary statistics from eqtlbma by Sarah:
	*	 [Midway] /project/mstephens/gtex/analysis/june2014/sumstats
*	 GTEx genotypes
	*	 /project/mstephens/gtex/preprocessing/june2014/omni451.vcf.gz
*	 GTEx expression level
	*	 [Midway] /project/mstephens/gtex/preprocessing/june2014/GTEx_Analysis_2014-06-13_expression_provisional/

## Release V6
Input and summary statistics from `MatrixEQTL`. See [this page](https://github.com/stephenslab/lab-resource/blob/master/data/gtex-v6-eqtl.md) for details.

### Why covariates not removed from expression values
"One reason for this is because this is how Matrix eQTL accepts input data. I realize that corrected expression values could be valuable to people and would certainly be willing to provide them, but I do not have a method. Although PEER does output a residuals file, I can only get PEER to run robustly by filtering low expressed genes and opting not to use explicit covariates as input to PEER. The solution would be to use the transformation that Matrix eQTL uses internally. ..." [David DeLuca, Broad institute]

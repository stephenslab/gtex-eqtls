# Computational Steps for Calculating Summary Statistics Under Quasi-Poisson Model
## First Attempt
Here is the command which breaks the data into 2,000 batches and perform analysis under quasi-Poisson model.

```bash
  nBatches=2000
  Model=quasipoisson
  for i in `seq $nBatches`; do
  echo '#!/bin/bash
       source $HOME/GIT/type-model/conf/GTEx.bashrc
       python $SrcDir/analysis_admin.py eqtlbma_batch \
       -g $InputDir/tss_coords.bed.gz \
       -s $InputDir/snp_coords.bed.gz \
       -n '"$nBatches"' -b '"$i"' -e ~/software/bin/eqtlbma_bf \
       -a $ConfDir/eqtlbma.'"$Model"'.txt' |\
  sbatch -J eqtlbma_$Model_$nBatches_$i -o $LogDir/eqtlbma_$Model_$nBatches_$i.o%j --mem-per-cpu=10000
  done
```

### Result
Computation is heavy and not completed (see next section). Here taking a sample batch out of the 2000 batches above, the following is summarized from `eqtlbma` log file

*	 Input batch has 28 genes and 113 SNPs
*	 Genotype and expression of 44 tissues are loaded, vary in sample size, from >70 to 450.
*	 Tissue specific genes "to keep" vary, but are around 26, 27, 28.
*	 Number of SNPs (from `omni451.vcf.gz`) drop greatly after discarding missing values and filter by MAF < 0.05. There are only 10 out of 113 SNPs left.
*	 Test for each gene-snp pair are done for every batch. There are only a handful such pairs for 2000 batches.

When 5000 batches are created in another test run, some batch may end up fail after loading the data because there will be no qualified SNP left for that batch, from our sample data.

### A note on computational time
There will be a 20min "overhead" for each job (loading / extracting data) regardless of per batch size, i.e., we will "waste" 0.3N CPU hours for running in N batches. It is thus best to find a batch neither too small nor to large (if we do not optimize the way data is loaded).

In my test job when only 10 SNPs and <28 genes are analyzed in effect, it takes only 1.38 sec to complete for normal model, evaluating for each tissue about 200 gene-SNP pairs.

However for the same job under quasi-Poission it takes **significantly** computation time. For example 5 hours (10:30 am to 3:30 pm) have passed for the 10 SNPs and 28 genes. Still the progress bar is at 3.75%. I have tried `--thread` option of `eqtlbma_bf` and reserve corresponding number of nodes per job on Midway. Somehow the computational node still only use one thread for a job. May have to figure it out but I doubt it's going to help now that we can easily create arbitrary batches.

Here is the `eqtlbma_bf` command in action for this one batch. For this batch `genes_2000_85_112.bed.gz` has 28 lines and `snps_2000_85_112.txt.gz` has 113 lines.

```bash
  /home/gaow/software/bin/eqtlbma_bf --geno /project/mstephens/gtex/analysis/june2014/types_v1/list_geno.txt \
                                     --exp /project/mstephens/gtex/analysis/june2014/types_v1/list_explevels.txt \
                                     --out eqtlbma_bf_quasipoisson_2000_85_112/eqtlbma_bf_quasipoisson_2000_85_112 \
                                     --anchor TSS --cis 100000 --lik quasipoisson --analys join --outss --maf 0.05 \
                                     --gridL /project/mstephens/gtex/analysis/june2014/types_v1/grid_phi2_oma2_general.txt.gz \
                                     --gridS /project/mstephens/gtex/analysis/june2014/types_v1/grid_phi2_oma2_with-configs.txt.gz \
                                     --bfs sin --error uvlr --fiterr 0.5 -v 1 \
                                     --gcoord genes_2000_85_112.bed.gz --snp snps_2000_85_112.txt.gz
```

## What Next?
After discussion with Matthew and William on June 22, 2015 we decide to pause this analysis (from experience of Pritchard's lab the gain from using GLM does not worth it's computational cost).

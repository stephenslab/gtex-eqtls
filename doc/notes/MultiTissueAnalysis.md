# Analyzing GTEx Using `eqtlbma`
## Pre-processing
### Extract input data from archive
```
snakemake extract_data
```

**_Important_**

 This is a large data-set. It is always not good idea to run any software / methods in development in this production-level data. My strategy for this project is to make a small subset of data to test on each step, fix the software / workflow and establish the procedure. Once established, I will run the analysis on the entire data-set. **Only comments / results fro the full data-set is documented here**.
 The subset of data is generated using the [same data](SumstatsDB.md) I have prepared for the `matrix-ash` method development. This also makes these two methods developments comparable.
 The `snakemake` workflow are identical between the subset and the full data, except the difference in `config.yaml` which is a symbolic link to configurations under either version. Switching between versions is just `make lite` vs `make full`.
 The Lite version is prepared by commands below, under `preprocessing` folder:
```bash
  snakemake get_lite_list
  snakemake prepare_lite_covariates
  snakemake prepare_lite_expression --config IsCluster=T
  snakemake prepare_lite_snps --config IsCluster=T
  cd /project/mstephens/gtex/analysis/april2015/Lite
  ln -s data-null data-dir
  ln -s data-max/GTEx_Analysis_2015-01-12_eQTLInputFiles_geneLevelNormalizedExpressionMatrices data-null/GTEx_Analysis_2015-01-12_eQTLInputFiles_geneLevelNormalizedExpressionMatrices
```

 I verified the data integrity: 16069 genes, 15507 unique max SNPs and 47643 unique null SNPs, matching with the lite data without an error. This is exactly the subset of data using which `matrix-ash` was developed.

### Generate Gene TSS / SNP coordinates files
The coordinate file for genes needs to be prepared. Also unfortunately the release does not come with a list of SNPs (union) involved. Getting such a list is quite a heavy duty. Takes 10min to extract ID's in parallel and 50min to concatenate them into a unique list in bed format.

```
snakemake prepare_coords
```

**_Note_**

 In Sarah's analysis based on an earlier version of data, there are 55,993 genes and 6,856,776 SNPs provided in gene/snp lists. In the v6 release there are 56,318 genes and 10,297,646 SNPs listed; although the sumstats of v6 only has ~39,000 genes.

### Generate input file lists
```
snakemake prepare_input_lists
```

## The Configuration Model
### Run `eqtlbma_bf` analysis in batches
`analysis_admin eqtlbma_batch` generates batches on the fly and perform `eqtlbma_bf` analysis in (embarrassing) parallel fashion on batches of genes. To use it,

```
python analysis_admin.py eqtlbma_batch -h
```

```
  usage: analysis_admin.py eqtlbma_batch [-h] -g GENE_COORDS -s SNP_COORDS
                                         [-w N] [-n N_BATCHES] [-b BATCH_ID] -a
                                         ARGS_FILE [--seed N] -e EQTLBMA_PATH
                                         [--dry-run] [--clean]
  
  optional arguments:
    -h, --help       show this help message and exit
    -g GENE_COORDS   Gene (TSS) coordinate file, in bed.gz format. (default:
                     None)
    -s SNP_COORDS    SNP coordinate file, in bed.gz format. (default: None)
    -w N             Window size. (default: 100000)
    -n N_BATCHES     Total number of batches. (default: 1000)
    -b BATCH_ID      Execute the i-th batch. Program will quit if invalid ID is
                     provided. (default: 1)
    -a ARGS_FILE     Path to file containing additional eqtlbma command
                     arguments. (default: None)
    --seed N         If specified, a random number will be generated using (N +
                     batch_id) as seed, and the eqtlbma command will be appended
                     a "--seed" argument with the number generated here.
                     (default: None)
    -e EQTLBMA_PATH  Path to an eqtlbma_* executable. (default: None)
    --dry-run        Only generate and save the batch data & commands without
                     performing analysis. (default: False)
    --clean          Remove batch gene / snp coords file upon finishing the
                     analysis. (default: False)
  
```

Here is the command which breaks the data into 100 batches and perform analysis

```bash
  nBatches=100
  Model=normal
  for i in `seq $nBatches`; do
  echo '#!/bin/bash
       source $HOME/GIT/gtex-eqtls/src/cfg/GTEx.bashrc
       python $SrcDir/python/analysis_admin.py eqtlbma_batch \
       -g $InputDir/tss_coords.bed.gz \
       -s $InputDir/snp_coords.bed.gz \
       -n '"$nBatches"' -b '"$i"' -w 100000 \
       --seed 10086 -e ~/software/bin/eqtlbma_bf \
       -a $ConfDir/eqtlbma.'"$Model"'.txt' |\
  sbatch -J eqtlbma_"$Model"_"$nBatches"_"$i" -e $LogDir/eqtlbma_"$Model"_"$nBatches"_"$i".e%j \
         -o $LogDir/eqtlbma_"$Model"_"$nBatches"_"$i".o%j --mem-per-cpu=10000 --time=36:00:00
  sleep 1
  done
```

See my lab notebook for various failures and fixes running this command. Read on for procedure that succeeded.

### Configuration model without permutation
Job submission

```bash
  nBatches=100
  Model=normal.perm0
  for i in `seq $nBatches`; do
  echo '#!/bin/bash
       source $HOME/GIT/gtex-eqtls/src/cfg/GTEx.bashrc
       python $SrcDir/python/analysis_admin.py eqtlbma_batch \
       -g $InputDir/tss_coords.bed.gz \
       -s $InputDir/snp_coords.bed.gz \
       -n '"$nBatches"' -b '"$i"' -w 100000 \
       --seed 10086 -e ~/software/bin/eqtlbma_bf \
       -a $ConfDir/eqtlbma.'"$Model"'.txt' |\
  sbatch -J eqtlbma_"$Model"_"$nBatches"_"$i" -e $LogDir/eqtlbma_"$Model"_"$nBatches"_"$i".e%j \
         -o $LogDir/eqtlbma_"$Model"_"$nBatches"_"$i".o%j --mem-per-cpu=10000 --time=36:00:00
  sleep 1
  done
```

Run completed under 24 hours, yielding to 45GB output file in various batches. For storage purpose it is best idea to first convert these files in different batches to HDF5 and merge them to one file & archive. Skipping this step for now as I want to move on to the HM and BMA steps. These output are archived as is, though:

```
tar -cvf eqtlbma_bf.normal.perm0.tar *
```

### List of summary statistics
The following script produces sumstats lists for all results previously computed:

```bash
  SumstatsDir=$ArchiveDir/eqtlbma_bf/July2015
  for i in `find $SumstatsDir -maxdepth 1 -name "eqtlbma_bf_normal_*" -type d`; do
      for j in `ls $i/*_sumstats_*.txt.gz`; do
          echo `echo $j | sed 's/\(.*sumstats_\)\(.*\).txt.gz/\2/g'` $j >> $i.ss.list
      done
  done
```

Future analysis can be based on summary statistics instead.

## The Type Model
The raw BF data are fed to the type-model version of the hierarchical model EM algorithm to estimate hyperparameters. Data is 28GB compressed text which may exceed 256GB when all are loaded into RAM (which is what `eqtlbma_hm` does!).

See my labnotes of July, 2015 to August, 2015 for failures and fixes.

# Analyzing GTEx Using `eqtlbma`
## Pre-processing
### Extract input data from archive
```bash
  for i in GTEx_Analysis_2015-01-12_eQTLInputFiles_covariates.tar.gz GTEx_Analysis_2015-01-12_eQTLInputFiles_snpMatrices.tar.gz GTEx_Analysis_2015-01-12_eQTLInputFiles_geneLevelNormalizedExpressionMatrices.tar.gz GTEx_Analysis_2015-01-12_eQTLInputFiles_snpMatricesSupplement.tar.gz; do
  for j in `tar -tf $i`; do
  echo $i $j
  mkdir -p $(basename $i .tar.gz)
  echo -e '#!/bin/bash\ntar -zxvf '"$i"' '"$j"' -O | gzip --best > '"$(basename $i .tar.gz)"'/'"$j"'.gz\necho complete!' | sbatch -J Extract."$i"."$j" -o $LogDir/Extract."$i"."$j".o%j
  done
  done
```

### Gene TSS coordinates file
```
GetTSSCoords "$DataPrefix"eQTLInputFiles_genePositions.txt.gz $InputDir/tss_coords.bed.gz
```

### SNP coordinates file
Unfortunately the release does not come with a list of SNPs (union) involved. Getting such a list is quite a heavy duty. Takes 10min to extract ID's in parallel and 50min to concatenate them into a unique list in bed format:

```
GetSNPCoords "$DataPrefix"eQTLInputFiles_snpMatrices $InputDir/snp_coords.bed.gz 1
GetSNPCoords "$DataPrefix"eQTLInputFiles_snpMatrices $InputDir/snp_coords.bed.gz 2
```

**_Note_**

 In Sarah's analysis based on an earlier version of data, there are 55,993 genes and 6,856,776 SNPs provided in gene/snp lists. In the v6 release there are 56,318 genes and 10,297,646 SNPs listed; although the sumstats of v6 only has ~39,000 genes.

### Input file lists
```bash
  rm -rf $InputDir/list_geno.txt
  DFolder="$DataPrefix"eQTLInputFiles_snpMatrices
  for i in `ls $DFolder`; do
      j=`basename $i _Analysis.snps.txt.gz`
      echo -e "$j\t$DFolder/$i" >> $InputDir/list_geno.txt
  done
  #
  rm -rf $InputDir/list_expr.txt
  DFolder="$DataPrefix"eQTLInputFiles_geneLevelNormalizedExpressionMatrices
  for i in `ls $DFolder`; do
      j=`basename $i _Analysis.expr.txt.gz`
      echo -e "$j\t$DFolder/$i" >> $InputDir/list_expr.txt
  done
  #
  rm -rf $InputDir/list_covar.txt
  DFolder="$DataPrefix"eQTLInputFiles_covariates
  for i in `ls $DFolder`; do
      j=`basename $i _Analysis.covariates.txt.gz`
      echo -e "$j\t$DFolder/$i" >> $InputDir/list_covar.txt
  done
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
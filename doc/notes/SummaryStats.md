# GTEx Summary Statistics
Computations are done on the Midway cluster. `GTEx.bashrc` sets up the proper shell environment for this computational workflow:

```
source $HOME/GIT/gtex-eqtls/src/cfg/GTEx.bashrc
```

## Formatting GTEx Summary Statistics V6
### Convert tissue specific results to matrices
For ease of storage / access / numeric operations I convert the summary statistics from plain text file to HDF5 format.

```bash
  cd $AnlyDir/GTEx_Analysis_2015-01-12_MatrixEQTL_allCisSNPGenePairs
  for i in `ls *.gz`; do
  echo '#!/bin/bash
       source $HOME/GIT/gtex-eqtls/src/cfg/GTEx.bashrc
       python $SrcDir/python/analysis_admin.py ss_to_h5 \
       '"$i"' --action convert --output . --message "MatrixEQTL summary statistics of GTEx Release 2015.04.15"' |\
  sbatch -J ss2h5_$i -o $LogDir/ss2h5_$i.o%j -e $LogDir/ss2h5_$i.e%j
  done
```

**_Note_**

 I ended up having to perform the entire data conversion on SSD on my desktop computer using. The `pytable` version of HDF5 implementation takes 2h to convert each tissue but on Midway when HDF5 file size grows to over 1GB the disk I/O is unacceptable. There must be tweaks but not looking into it for the moment. Implementations with `h5py` for the same purpose is way faster in data conversion (30min each file, performance on Midway is also acceptable) but the resulting data file size is 1.4 times of `pytables` implementation, even both set to zlib9 compression -- thus I stick with `pytable`. Once converted, it takes almost zero seconds to query data from the HDF5 files.

### Merge data across tissues
```bash
  source $HOME/GIT/gtex-eqtls/src/cfg/GTEx.bashrc
  zcat $InputDir/tss_coords.bed.gz | cut -f4 | split -l 1127 -d - $TmpDir/tss_coords.bed.split.
  k=0
  for i in `ls $TmpDir/tss_coords.bed.split.*`; do
  let k+=1
  echo '#!/bin/bash
       source $HOME/GIT/gtex-eqtls/src/cfg/GTEx.bashrc
       python $SrcDir/python/analysis_admin.py ss_to_h5 \
       $DataDir/GTEx_Analysis_2015-01-12_MatrixEQTL_allCisSNPGenePairs/*.h5 \
       --action merge --output $DataDir/MatrixEQTLSumStats/MatrixEQTLSumStats.'"$k"' \
       --message "MatrixEQTL summary statistics of GTEx Release 2015.04.15" \
       --gene-list '"$i"'' |\
  sbatch -J h5merge_"$k" -o $LogDir/h5merge_"$k".o%j -e $LogDir/h5merge_"$k".e%j
  done
```

It takes on average 3 hours per batch to complete. Merged number of genes per batch see `log/2015-06-18-merger-log.txt`. There are a total of 38933 genes.

**_Note_**

 I was suspicious to see only 38933 genes in v6. I decide to check against the input data as well as the intermediate HDF5 files I generated
```bash
  for i in `ls $DataDir/GTEx_Analysis_2015-01-12_MatrixEQTL_allCisSNPGenePairs/*.cis.eqtl.gz`; do
  echo '#!/bin/bash
       source $HOME/GIT/gtex-eqtls/src/cfg/GTEx.bashrc
       zcat '"$i"' | cut -f2 | sort -u | gzip --best > '"$i"'.uniqgenes.gz' |\
  sbatch -J count_gene -o $LogDir/count_gene.o%j -e $LogDir/count_gene.e%j
  done
  for i in `ls $DataDir/GTEx_Analysis_2015-01-12_MatrixEQTL_allCisSNPGenePairs/*.h5`; do
  echo '#!/bin/bash
       source $HOME/GIT/gtex-eqtls/src/cfg/GTEx.bashrc
       python $SrcDir/python/analysis_admin.py ss_to_h5 '"$i"' --action summary --output stdout | gzip --best > '"$i"'.h5genes.gz' |\
  sbatch -J count_gene_h5 -o $LogDir/count_gene_h5.o%j -e $LogDir/count_gene_h5.e%j
  done
```

 Counting (`zcat *.xx.gz | sort -u | wc -l`) and comparing results from both sources, I found there are indeed only 38933 genes in the output. So we are good.

Now merging the batches into one large HDF5 file

```
python $SrcDir/python/analysis_admin.py ss_to_h5 $DataDir/MatrixEQTLSumStats/*.h5 --action cat --output $DataDir/MatrixEQTLSumStats.h5
```

### Select the "best" gene-snp pair
```
python $SrcDir/python/analysis_admin.py ss_to_h5 $DataDir/MatrixEQTLSumStats.h5 --action max --output $DataDir/MatrixEQTLSumStats.Portable.h5
```

### Select some "null" gene-snp pairs
For each gene I select up to 3 null gene-snp pairs, for use of calibrating Sarah's model

```
python $SrcDir/python/analysis_admin.py ss_to_h5 $DataDir/MatrixEQTLSumStats.h5 --action null --output $DataDir/MatrixEQTLSumStats.Portable.h5 --gene-list <(zcat $DBDir/snp-gene-pairs/GTExV6Genes.gz)
```

## Creating R Interface for the New Format
Since the release did not use rsID to name SNPs I create a database to match SNP ID with dbSNP names, to use with R to extract proper SNP-gene by rsID.

```
wget ftp://ftp.ncbi.nih.gov/snp/organisms/human_9606_b144_GRCh37p13/VCF/All_20150605.vcf.gz
```

Build an annotation database matching SNP ID to rs ID as well as the cis-genes involved of up to \(\pm100,000\)bp.

```bash
  source $HOME/GIT/gtex-eqtls/src/cfg/GTEx.bashrc
  zcat $InputDir/snp_coords.bed.gz | split -l 205953 -d - $TmpDir/snp_coords.bed.split.
  k=0
  for i in `ls $TmpDir/snp_coords.bed.split.*`; do
  let k+=1
  echo '#!/bin/bash
       source $HOME/GIT/gtex-eqtls/src/cfg/GTEx.bashrc
       cat '"$i"' | bash $SrcDir/shell/FindrsIDcisGene.sh \
       $DBDir/dbSNP/All_20150605.vcf.gz $InputDir/tss_coords.bed.gz 100000 | gzip > \
       $DBDir/snp-gene-pairs.'"$k"'.gz' |\
  sbatch -J snpid_"$k" -o $LogDir/snpid_"$k".o%j -e $LogDir/snpid_"$k".e%j
  done
```

It takes on average 25 hours per batch to complete. There are 9794339 out of 10297646 (95.1%) SNPs with rsID in latest dbSNP (build 144). To merge them to one file and build a sqlite database:

```
zcat $DBDir/snp-gene-pairs.*.gz | gzip --best > $DBDir/snp-gene-pairs.gz
cd $DBDir; mkdir snp-gene-pairs; mv snp-gene-pairs.* snp-gene-pairs; cd snp-gene-pairs
wsqlite snp-gene.sqlite3 -i snp-gene-pairs.gz --as dbsnp144 -d '\t' --header coord rsid cisgenes
wsqlite snp-gene.sqlite3 "create index rsid_index on dbsnp144 (rsid)"
wsqlite snp-gene.sqlite3 "create index coord_index on dbsnp144 (coord)"
```

### Run queries
Now we have the data `*.h5` and meta-data `*.sqlite` in place. I have made a separate note on how to use this data-set. See `SumstatsDB.pdf` in this repo.

## Reproducing GTEx V6 Cis-eQTL Analysis
The idea now is to figure out how the `matrixEQTL` analysis is done exactly in v6 so that `eqtlbma` will be done on the same basis, i.e., input files and parameters. For this purpose I make a toy dataset for a randomly chosen tissue from GTEx consisting of five files for input to `matrixEQTL`: genotype `SNP.txt`, expression `GE.txt`, a file `Covariates.txt` of covariates and files `geneloc.txt` and `snpsloc.txt` with gene and SNP location information.

```
Tissue=Adipose_Subcutaneous
```

### Covariates
There are 40 covariates: "C1,C2,C3,InferredCov1, InferredCov2, ... InferredCov35, gender, Platform". Not all are included in every analysis -- inclusion depends on sample size. See the `README` for details.

```
tar -zxvf GTEx_Analysis_2015-01-12_eQTLInputFiles_covariates.tar.gz "$Tissue"_Analysis.covariates.txt
wc -l "$Tissue"_Analysis.covariates.txt
cut -f1 "$Tissue"_Analysis.covariates.txt | tr "\n" ", "
mv "$Tissue"_Analysis.covariates.txt $TmpDir/Covariates.txt
```

### SNP data
#### SNP matrix
Input genotype data are large. Here I just take the first 500 SNPs from the tissue data for verification purpose.

```
tar -zxvf GTEx_Analysis_2015-01-12_eQTLInputFiles_snpMatrices.tar.gz "$Tissue"_Analysis.snps.txt -O |\
    head -501 > $TmpDir/SNP.txt
```

**_Warning_**

 SNP data released has a lot of imputed data. Cautions are to be taken when comparing with previous `eqtlbma` analysis (have to make sure the input data are the same).

#### SNP location
This is required input for `matrixEQTL` but is not provided. Need to create it from SNP ID's.

```
echo -e "snp\tchr\tpos" > $TmpDir/snpsloc.txt
cut $TmpDir/SNP.txt -f1 | tail -n+2 | awk -F"_" '{print $0"\t""chr"$1"\t"$2}' >> $TmpDir/snpsloc.txt
```

### Expression data
To focus only on genes paired with the 500 SNPs I get a list of genes via `bedtools closest` and only extract expression data for these genes.

Create `bedtools` input and run `bedtools` to search for genes within 1000bp range to the SNPs of interest.

```
tail -n+2 GTEx_Analysis_2015-01-12_eQTLInputFiles_genePositions.txt |\
     awk '{print $2"\t"$3"\t"$4+1"\t"$1"\t"1"\t""+"}' > $TmpDir/geneloc.bed
cut $TmpDir/SNP.txt -f1 | tail -n+2 | awk -F"_" '{print $1"\t"$2"\t"$2+1"\t"$0"\t"1"\t""+"}' > $TmpDir/snpsloc.bed
bedtools closest -a $TmpDir/snpsloc.bed -b $TmpDir/geneloc.bed |\
         awk '{if (($8-$2)<1000 && ($8-$2)>-1000) print $0}' > $TmpDir/gs-pairs.bed
```

Extract expression data for selected genes

```
cut -f10 $TmpDir/gs-pairs.bed | sort -u > $TmpDir/1.txt # 22 genes found
tar -zxvf GTEx_Analysis_2015-01-12_eQTLInputFiles_geneLevelNormalizedExpressionMatrices.tar.gz \
    "$Tissue"_Analysis.expr.txt -O | head -1 > $TmpDir/GE.txt
tar -zxvf GTEx_Analysis_2015-01-12_eQTLInputFiles_geneLevelNormalizedExpressionMatrices.tar.gz \
    "$Tissue"_Analysis.expr.txt -O | grep -wf $TmpDir/1.txt >> $TmpDir/GE.txt # 11 genes here
```

It ends up with 11 genes for the 500 SNPs

```
cut -f 1,2 $TmpDir/GE.txt
```

### Gene location data
Location data is extracted for the selected genes.

```
cut -f1 $TmpDir/GE.txt | tail -n+2 > $TmpDir/1.txt
grep -wf $TmpDir/1.txt GTEx_Analysis_2015-01-12_eQTLInputFiles_genePositions.txt > $TmpDir/2.txt
head -1 GTEx_Analysis_2015-01-12_eQTLInputFiles_genePositions.txt | cat - $TmpDir/2.txt > $TmpDir/geneloc.txt
```

### Find overlapping gene-snp pairs between this toy set and GTEx v6 summary statistics
```bash
  for i in `cut -f1 $TmpDir/geneloc.txt | tail -n+2`; do
      echo $i
      python analysis_admin.py sumstat_query $InputDir/"$Tissue"_Analysis.h5 -g $i -s NULL > $TmpDir/1.txt
      if [[ -s $TmpDir/1.txt ]]; then
         grep -wf $TmpDir/1.txt $TmpDir/snpsloc.txt | cut -f1 > $TmpDir/"$i"_snps_to_check.txt
      fi
  done
```

### Run `matrixEQTL` on this toy set
Once the input is prepared it is straightforward to just download the example R script (below, also committed to git repo) and run it, just remember to set the correct data path and change emitted p-value threshold to 1 (`pvOutputThreshold=1`).

http://www.bios.unc.edu/research/genomic_software/Matrix_eQTL/R.html

To obtain value for one gene-snp pair:

```r
  source('MatrixEQTL.R')
  ss = '1_714019_A_G_b37'
  gg = 'ENSG00000237491.4'
  me$all$eqtls[which(me$all$eqtls$snps == ss & me$all$eqtls$gene == gg ),]
```

```
                  snps              gene statistic       pvalue          FDR
  241 1_714019_A_G_b37 ENSG00000237491.4  5.392883 1.576441e-07 3.597687e-06
           beta
  241 0.7902205
```

and to verify value for gene-snp pairs:

```
echo $Tissue
```

```r
  source('/project/mstephens/gtex/scripts/SumstatQuery.R')
  dat <- GetFlatSS('ENSG00000237491.4', '/project/mstephens/gtex/analysis/april2015/eqtl_data/GTEx_Analysis_2015-01-12_MatrixEQTL_allCisSNPGenePairs/Adipose_Subcutaneous_Analysis.h5')
  print(dat["1_714019_A_G_b37", ])
```

```
          beta       t-stat      p-value
  7.902205e-01 5.392883e+00 1.576441e-07
```

I verified a few other gene-SNP pairs, all agree. So I can claim I've figured out the input data and parameters used for the GTEx analysis!

## Update `eqtlbma` Analysis for Some Gene-snp Pairs and Verify with GTEx `matrixEQTL`
Here is the command for generating summary statistics via `eqtlbma_bf`:

```bash
  Analysis=join
  Tissue=Test_Tissue
  echo -e "$Tissue\tSNP.txt.gz" > list_geno.txt
  echo -e "$Tissue\tGE.txt" > list_expl.txt
  echo -e "$Tissue\tCovariates.txt" > list_covar.txt
  eqtlbma_bf --geno list_geno.txt --exp list_expl.txt --covar list_covar.txt \
             --scoord snpsloc.bed  --gcoord geneloc.bed \
             --out ToyEQTLBma --analys $Analysis --outss --bfs sin \
             --gridL /project/mstephens/gtex/analysis/june2014/types_v1/grid_phi2_oma2_general.txt.gz \
             --gridS /project/mstephens/gtex/analysis/june2014/types_v1/grid_phi2_oma2_with-configs.txt.gz \
             --error uvlr --thread 4
```

Check a couple of output from `eqtlbma_bf` with `matrixEQTL`

```
zcat ToyEQTLBma_sumstats_Test_Tissue.txt.gz | head | cut -f1,2,7,9
```

```
  gene    snp     betahat.geno    betapval.geno
  ENSG00000177757.1       1_662622_G_A_b37        -2.588170e-04   9.985439e-01
  ENSG00000177757.1       1_676127_C_T_b37        7.151412e-02    7.744839e-01
  ENSG00000177757.1       1_691541_AT_A_b37       -9.077870e-02   4.383949e-01
  ENSG00000177757.1       1_693625_T_C_b37        -5.757815e-02   7.415218e-01
  ENSG00000177757.1       1_693731_A_G_b37        -1.080375e-02   9.255429e-01
  ENSG00000177757.1       1_697411_G_GA_b37       5.502301e-01    2.327584e-02
  ENSG00000177757.1       1_701835_T_C_b37        7.243649e-01    3.168885e-03
  ENSG00000177757.1       1_704367_T_C_b37        1.898470e-01    3.756234e-01
  ENSG00000177757.1       1_705882_G_A_b37        -1.595642e-01   3.141656e-01
```

Check with `MatrixEQTL`

```r
  gg = 'ENSG00000177757.1'
  for (ss in c('1_662622_G_A_b37', '1_676127_C_T_b37', '1_691541_AT_A_b37', '1_693625_T_C_b37', '1_693731_A_G_b37')) print(me$all$eqtls[which(me$all$eqtls$snps == ss & me$all$eqtls$gene == gg ),][,c(2,1,6,4)])
```

```
                    gene             snps         beta    pvalue
  5493 ENSG00000177757.1 1_662622_G_A_b37 -0.000258817 0.9985439
                    gene             snps       beta    pvalue
  4593 ENSG00000177757.1 1_676127_C_T_b37 0.07151412 0.7744839
                    gene              snps       beta    pvalue
  3047 ENSG00000177757.1 1_691541_AT_A_b37 -0.0907787 0.4383949
                    gene             snps        beta    pvalue
  4454 ENSG00000177757.1 1_693625_T_C_b37 -0.05757815 0.7415218
                    gene             snps        beta    pvalue
  5200 ENSG00000177757.1 1_693731_A_G_b37 -0.01080375 0.9255429
```

Thus `eqtlbma_bf` agrees with `MatrixEQTL` when input agree. Good! We'll have to update the analysis using new input data.
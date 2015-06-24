source('/project/mstephens/gtex/scripts/SumstatQuery.R')
##
# Look up GTEx SNP ID
##
ShowSNP('rs6600419', '/project/mstephens/gtex/analysis/april2015/query/snp-gene.db')
## GTEx SNP ID: 1_43124701_A_G_b37 
## cisGenes: ENSG00000065978.13,ENSG00000164007.6,ENSG00000171960.6,ENSG00000200254.1,ENSG00000234917.1,ENSG00000236180.2
##
# Load data for given gene
##
dat <- GetSS('ENSG00000171960.6', '/project/mstephens/gtex/analysis/april2015/query/MatrixEQTLSumStats.h5')
print(names(dat))
## [1] "beta"    "p-value" "t-stat"  "z-score"
##
# Output information for the SNP of interest
##
print(dat$"p-value"["1_43124701_A_G_b37", ])
 ##                 Adipose_Subcutaneous              Adipose_Visceral_Omentum 
 ##                         1.632499e-01                          3.462833e-02 
 ##                        Adrenal_Gland                          Artery_Aorta 
 ##                         8.168959e-01                          5.600338e-01 
 ##                      Artery_Coronary                         Artery_Tibial 
 ##                         7.659426e-01                          7.397466e-01 
 ## Brain_Anterior_cingulate_cortex_BA24           Brain_Caudate_basal_ganglia 
 ## ...
print(dat$"z-score"["1_43124701_A_G_b37", ])
 ##                 Adipose_Subcutaneous              Adipose_Visceral_Omentum 
 ##                          -1.39422418                           -2.11267803 
 ##                        Adrenal_Gland                          Artery_Aorta 
 ##                          -0.23153929                            0.58279135 
 ##                      Artery_Coronary                         Artery_Tibial 
 ##                           0.29768626                            0.33218889 
 ## Brain_Anterior_cingulate_cortex_BA24           Brain_Caudate_basal_ganglia
 ## ... 
###
# Load the best gene-snp data
###
dat <- GetSS('max', '/project/mstephens/gtex/analysis/april2015/query/MatrixEQTLSumStats.Max.h5')
print(dat$"p-value"[1, ])
print(dat$"z-score"[1, ])

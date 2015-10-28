source("/project/mstephens/gtex/scripts/SumstatQuery.R")
##
# Load data for given gene
##
dat <- GetSS("ENSG00000171960.6", "/project/mstephens/gtex/analysis/april2015/query/MatrixEQTLSumStats.h5")
print(names(dat))
## [1] "beta"    "p-value" "t-stat"  "z-score"
dim(dat$"beta")
## [1] 6344   44
dat$"beta"[1:4,1:4]
##                     Adipose_Subcutaneous Adipose_Visceral_Omentum Adrenal_Gland Artery_Aorta
## 1_42124161_C_T_b37            0.04608095                      NaN           NaN   0.06785597
## 1_42124246_A_AC_b37          -0.03087348              -0.03222961    -0.1286249  -0.13490575
## 1_42124311_G_A_b37            0.12400229               0.08011222           NaN  -0.21858459
## 1_42124614_C_T_b37            0.04630340                      NaN           NaN   0.06785597
##
# Output information for the SNP of interest
##
dat$"p-value"["1_43124701_A_G_b37", ]
 ##                 Adipose_Subcutaneous              Adipose_Visceral_Omentum 
 ##                         1.632499e-01                          3.462833e-02 
 ##                        Adrenal_Gland                          Artery_Aorta 
 ##                         8.168959e-01                          5.600338e-01 
 ##                      Artery_Coronary                         Artery_Tibial 
 ##                         7.659426e-01                          7.397466e-01 
 ## Brain_Anterior_cingulate_cortex_BA24           Brain_Caudate_basal_ganglia 
 ## ...
dat$"z-score"["1_43124701_A_G_b37", ]
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
mdat <- GetSS("max", "/project/mstephens/gtex/analysis/april2015/query/MatrixEQTLSumStats.h5")
dim(mdat$"t-stat")
mdat$"p-value"[1:4,1:4]
mdat$"t-stat"["ENSG00000000419.8_20_49461813_G_C_b37",]
## Is this gene-snp pair most significant in spleen?
dat <- GetSS("ENSG00000000419.8", "/project/mstephens/gtex/analysis/april2015/query/MatrixEQTLSumStats.h5")
idx.to.show <- matxMax(abs(dat$"t-stat"))
rownames(dat$"t-stat")[idx.to.show[1]]
colnames(dat$"t-stat")[idx.to.show[2]]
###
# Load the "null" gene-snp data
###
ndat <- GetSS("null", "/project/mstephens/gtex/analysis/april2015/query/MatrixEQTLSumStats.h5")
dim(ndat$"t-stat")
##
# Look up GTEx SNP ID
##
ShowSNP("rs6600419", "/project/mstephens/gtex/analysis/april2015/query/snp-gene.db")
## GTEx SNP ID: 1_43124701_A_G_b37 
## cisGenes: ENSG00000065978.13,ENSG00000164007.6,ENSG00000171960.6,ENSG00000200254.1,ENSG00000234917.1,ENSG00000236180.2
##
# Look up rs ID given GTEx SNP ID
##
ShowSNP("1_43124701_A_G_b37", "/project/mstephens/gtex/analysis/april2015/query/snp-gene.db")
## rsID(s): rs6600419 
## cisGenes: ENSG00000065978.13,ENSG00000164007.6,ENSG00000171960.6,ENSG00000200254.1,ENSG00000234917.1,ENSG00000236180.2 
##
# Create matched training/testing sets
##
N1 <- 8000
N2 <- 16069
strong.train <- SubsetMatLists(mdat, seq(1, N1))
strong.test <- SubsetMatLists(mdat, seq(N1 + 1, N2))
strong.train.genes <- as.character(lapply(strsplit(rownames(strong.train$beta), "_"), function(x) x[1]))
strong.test.genes <- as.character(lapply(strsplit(rownames(strong.test$beta), "_"), function(x) x[1]))
null.genes <- as.character(lapply(strsplit(rownames(ndat$beta), "_"), function(x) x[1]))
null.train <- SubsetMatLists(ndat, which(null.genes %in% strong.train.genes))
null.test <- SubsetMatLists(ndat, which(null.genes %in% strong.test.genes))

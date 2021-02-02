#!/usr/bin/Rscript

args <- commandArgs(TRUE);
BETAN <- args[1]
SEN <- args[2]
BETA <- args[3]
SE <- args[4]
OUT <- args[5]

source("/gpfs/commons/groups/lappalainen_lab/skim/fastqtl/functions_eQTLresult_II.R")
suppressPackageStartupMessages(library(mashr))

## read data
beta.n=read.table(BETAN, head=T, as.is=T, row.names=1); 
cat.dim(beta.n)
se.n=read.table(SEN, head=T, as.is=T, row.names=1); 
cat.dim(se.n)
beta=read.table(BETA, head=T, as.is=T, row.names=1); 
cat.dim(beta)
se=read.table(SE, head=T, as.is=T, row.names=1); 
cat.dim(se)

## estimate null correlation
data.temp = mash_set_data(beta.n, se.n, alpha=1)
Vhat = estimate_null_correlation(data.temp,  apply_lower_bound = FALSE)
rm(data.temp)

## make mash files
data.random = mash_set_data(beta.n, se.n, V=Vhat, alpha=1)
data.strong = mash_set_data(beta, se, V=Vhat, alpha=1)

## get data driven covariances
npc=ifelse(ncol(beta) < 5, ncol(beta), 5)
U.pca = cov_pca(data.strong,npc)
U.ed = cov_ed(data.strong, U.pca)

## fit the model
U.c = cov_canonical(data.random)
m = mash(data.random, Ulist = c(U.ed,U.c), outputlevel = 1)

## compute posterior summaries
m2 = mash(data.strong, g=get_fitted_g(m), fixg=TRUE)
print(length(get_significant_results(m2)));
print(length(get_significant_results(m2, thres=0.01)));
print(get_loglik(m2));
save(m2, file=paste0(OUT, "_m2.Rda"))

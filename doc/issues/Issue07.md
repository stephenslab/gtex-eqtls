## Issue 07: SVD Error of Eqtlbma_bf on GTEx V6 Data
This error occurs in 5% of tests on original data and 100% tests on permuted data; it prevents the program from going further:

```
  gsl: svd.c:286: ERROR: svd of MxN matrix, M<N, is not implemented
  Default GSL error handler invoked.
```

Need to identify the cause (create a toy data set) and fix the code

## Plan
For starters, the most difficult of this debug game is to find the particular gene-snp pair that can reproduce the issue. It is perhaps best to 1) run without permutation and monitor the progress bar to determine roughly which genes might be in question, and 2) extract all relevant data for those genes to create a smaller set, and 3) repeat 1) and 2) until identified the cause with a particular gene-snp pair data.

**_Note_**

 While trying to reproduce the issue on a toy data set I found another (or the more fundamental) bug of segfault that can be reproduced on `data/Issue07A.tar.gz`. The bug was fixed https://github.com/gaow/eqtlbma/commit/2c080acc2ea19c2070d67785caeaa53eb8a36aca.

## Resolution
A toy data was created `data/Issue07B.tar.gz` that reproduces the problem. The solution was to check if the input is proper before feeding it to GSL, see https://github.com/gaow/eqtlbma/commit/f17b2a5d34b80cc0be8fb53190c5818bde712388.
## Convergence Problem with the EM Algorithm for Type Model
### Problem
Fitting type model with 5% of data takes long time (>300 hours on 16 threads) yet eventually ended up not converging. EM algorithm halted because of decrease in likelihood.

### Diagnosis
```
for i in `ls -rt 2015{09*,10*}/*.log`; do grep loglik $i | awk '{print $1,$2,$4,$6}' | awk 'NR==1; END{print}' ; done
```

```
  iter 0 990600.446390 5.0000e-01
  iter 20 1042136.374448 2.7283e-04
  iter 0 1042136.526499 2.7283e-04
  iter 20 1044912.298293 1.2421e-06
  iter 0 1044912.189032 1.2421e-06
  iter 20 1045557.401620 1.9939e-09
  iter 0 1045557.552082 1.9939e-09
  iter 20 1045765.102952 2.3444e-12
  iter 0 1045765.253749 2.3444e-12
  iter 35 1045877.629320 1.5105e-17
  iter 0 1045877.753775 1.5105e-17
  iter 26 1045904.075885 2.0291e-21
  iter 0 1045904.129649 2.0291e-21
  iter 19 1045913.622268 2.9813e-24
  iter 0 1045917.080730 6.8150e-26
  iter 11 1045920.000938 1.5562e-27
  iter 0 1045919.943270 1.5562e-27
  iter 16 1045922.895941 6.3673e-30
  iter 0 1045922.803477 6.3673e-30
  iter 13 1045924.612012 7.2976e-32
  iter 0 1045924.798100 7.2976e-32
  iter 1 1045924.729358 5.1746e-32
```

Notice the last iteration from previous run has the same parameter as the initial iteration from the following run (in this text just the last column \(\pi_0\) are the same but in fact other parameters are also identical), which is expected. However the initial likelihood estimate differ. This is not right and should be purely numeric issue.
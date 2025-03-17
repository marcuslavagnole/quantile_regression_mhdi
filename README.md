In this repository, one can find the R routines used in the article [Bayesian variable selection in quantile regression with random effects: An application to Municipal Human Development Index](https://doi.org/10.1080/02664763.2021.1950654). It is a joint work with Kelly C. M. Gonçalves, published in the _Journal of Applied Statistics_.



This repo includes:

- **MCMC_SS-RE-GAL.R** : MCMC routine for the quantile regression with random effects and variable selection based on Generalized Asymmetric Laplace distribution;
- **MCMC_SS-RE-AL.R** : MCMC routine for the quantile regression with random effects and variable selection based on Asymmetric Laplace distribution;
- **MCMC_SS-RE-T.R** : MCMC routine for the linear regression with random effects and variable selection based on T distribution;
- **MCMC_SS-RE-N.R** : MCMC routine for the linear regression with random effects and variable selection based on Normal distribution;
- Directory **./Full_Conditionals** contains the full conditional distributions for the above models;
- **data.RData** : real data set.

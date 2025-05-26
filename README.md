This repository provides the R routines from the article [Bayesian variable selection in quantile regression with random effects: An application to Municipal Human Development Index](https://doi.org/10.1080/02664763.2021.1950654). It is a joint work with Kelly C. M. Gonçalves, published in the _Journal of Applied Statistics_.

According to the Atlas of Human Development in Brazil, the income dimension of the Municipal Human Development Index (MHDI-I) is an indicator that shows the population's ability in a municipality to ensure a minimum standard of living to provide their basic needs, such as water, food, and shelter. In public policy, one of the research objectives is to identify social and economic variables associated with this index. Due to the income inequality, evaluating these associations in quantiles instead of the mean could be more interesting. Motivated by the analysis of MHDI-I from municipalities in Rio de Janeiro, the paper develops a Bayesian variable selection in quantile regression models with hierarchical random effects. In particular, we assume a likelihood function based on the Generalized Asymmetric Laplace distribution, and a spike-and-slab prior is used to perform variable selection. The Generalized Asymmetric Laplace distribution is a more general alternative than the Asymmetric Laplace one, a common approach used in quantile regression under the Bayesian paradigm. 

The repo includes:

- **MCMC_SS-RE-GAL.R** : MCMC routine for the quantile regression with random effects and variable selection based on the Generalized Asymmetric Laplace distribution;
- **MCMC_SS-RE-AL.R** : MCMC routine for the quantile regression with random effects and variable selection based on the Asymmetric Laplace distribution;
- **MCMC_SS-RE-T.R** : MCMC routine for the linear regression with random effects and variable selection based on the T distribution;
- **MCMC_SS-RE-N.R** : MCMC routine for the linear regression with random effects and variable selection based on the Normal distribution;
- Directory **./Full_Conditionals** contains the full conditional distributions for the above models;
- Directory **./Auxiliary_Functions** contains auxiliary functions for running the above models;
- **data.RData** : real data set.

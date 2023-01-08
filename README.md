# R code for the Point Mass Proposal Metropolis-Hastings (PMPMH) algorithm

This repository contains code for general implementations and specific use cases of the PMPMH algorithm in [1]. The PMPMH algorithm uses a hidden Markov model (HMM) approximation to a state-space model (SSM) to provide an efficient proposal distribution for the latent states, allowing Metropolis-Hastings (M-H) sampling of the latent state vector and model parameters. 

## Files 
*PMPMH/functions/gen_blocks.R*: generates correctly formatted block indices for block updating of the latent states.

*PMPMH/functions/gen_HMM.R*: for a specified block of states, calculates the discrete HMM approximation to the SSM for a given set of tuning parameters.

*PMPMH/functions/gen_proposal.R*: samples a sequence of intervals from the HMM approximation and continuous values for the state at each time point from within those intervals.

*PMPMH/functions/run_MH.R*: runs one iteration of the M-H algorithm's accept/reject step.
 

*PMPMH/implementations/PMPMH_GaussMix_mod1_approach\*.R*: use case of PMPMH method for fitting Gaussian mixture SSM Model 1 (Model 1 of the first numerical illustration in [1]) under approach \* for defining the grid cells.

*PMPMH/implementations/PMPMH_GaussMix_mod2_approach\*.R*: use case of PMPMH method for fitting Gaussian mixture SSM Model 2 (Model 2 of the first numerical illustration in [1]) under approach \* for defining the grid cells.

*PMPMH/implementation/PMPMH_Nicholsons_approach\*.R*: use case of the PMPMH method for SSM fitting to data simulated from Nicholson's blowfly model (2-dimensional numerical illustration of [1]) under approach \* for defining the grid cells.
 

*PG/\*.R*: implementations of the particle Gibbs algorithm (both case studies in [1]). Code based on github.com/nchopin/particles.

*PGAS/\*.R*: implementations of the particle Gibbs with Ancestor sampling algorithm (both case studies in [1]). Code based on github.com/nchopin/particles and [2]. 


## Configuration/operation 
Clone repository to preferred installation of R.The PMPMH functions for general one-dimensional SSMs, and two-dimensional SSMs with conditionally updated states, are already included. For the implementation of examples not given in cloned implementation files, follow the implementation files as a template for supplied arguments. 

Required packages: (CRAN) `invgamma`, `truncnorm`. 
CRAN installation:

``` 
#install.packages("invgamma")
library(invgamma)
#install.packages("truncnorm")
library(truncnorm)
```

## Development
This code is maintained by Mary Llewellyn (mary.llewellyn@ed.ac.uk)

For troubleshooting or comments, please submit an issue or e-mail.

## Attributions
If you find the PMPMH algorithm useful for your research, you can acknowledge by citing [1].

## Update information

## References 
[1] 

[2] Lindsten, F., Jordan, M. I., and Schön, T. B. (2014). Particle Gibbs with Ancestor Sampling. Journal of Machine Learning Research, 15(63):2145–2184.



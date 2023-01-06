# R code for the Point Mass Proposal Metropolis-Hastings (PMPMH) algorithm

This repository contains code for general implementations and specific use cases of the PMPMH algorithm in [1]. The PMPMH algorithm uses a hidden Markov model (HMM) approximation to a state-space model (SSM) to provide an efficient proposal distribution for the latent states, allowing Metropolis-Hastings (M-H) sampling of the latent state vector and model parameters. 

## Files 
*functions/gen_blocks.R*: generates correctly formatted block indices for block updating of the latent states.

*functions/gen_HMM.R*: for a specified block of states, calculates the discrete HMM approximation to the SSM for a given set of tuning parameters.

*functions/gen_proposal.R*: samples a sequence of intervals from the HMM approximation and continuous values for the state at each time point from within those intervals.

*functions/run_MH.R*: runs one iteration of the M-H algorithm's accept/reject step.

*implementations/Mixture_of_Gaussians.R*: use case of PMPMH method for fitting Gaussian mixture SSM to data (case study 1 in [1]).

*implementation/Nicholsons.R*: use case of the PMPMH method for SSM fitting to data simulated from Nicholson's blowfly model (case study 2 in [1]).

## Configuration/operation 
Clone repository to preferred installation of R. Functions for general one-dimensional SSMs are already included. For implementation of examples not given in cloned implementation files, follow the implementation files as a template for supplied arguments. 

required packages: (CRAN) `invgamma`, `truncnorm`. CRAN installation:

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
[no updates: 21/03/2022]

## References 
[1] Llewellyn, M., King, R., Elvira, V., and Ross, G. (2022). A Point Mass Proposal Method for Bayesian State-Space Model Fitting. *https://doi.org/10.48550/arXiv.2203.13649* (Author's Original Manuscript)



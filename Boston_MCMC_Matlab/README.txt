Spatial econometrics models fit by R. S. Bivand.

The raw SE Matlab results for the Boston data set with 490 observations 
(the uncensored ones) and a queen weights matrix.

All 25000 draws - 5000 omitted,

sem_g is an mcmc object, 1:14 beta, 15 rho_{Err}, 16 sigma
sdem_g is an mcmc object, 1:27 beta, 28 rho_{Err}, 29 sigma
slx_g is an mcmc object, 1:27 beta, 28 sigma

sar_g is an mcmc object, 1:14 beta, 15 rho_{Lag}, 16 sigma, the file 
includes direct and total 20000x13 matrixes

sdm_g is an mcmc object, 1:27 beta, 28 rho_{Lag}, 29 sigma, the file 
includes direct and total 20000x13 matrixes

sac_g is an mcmc object, 1:14 beta, 15 rho_{Lag}, 16 rho_{Err}, 17 sigma, 
the file includes direct and total 20000x13 matrixes

sdac_g failed to run (Manski model).

The direct and total matrices could be put through as.mcmc() for easier 
summaries.

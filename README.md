# Spatial lag model with INLA: the `slm` latent effect

This repository contains the `R` code of paper [*Estimating Spatial Econometrics Models with Integrated Nested Laplace Approximation*](https://arxiv.org/abs/1703.01273) by Virgilio Gómez-Rubio, Roger S. Bivand, Håvard Rue.


Files available are:

* `boston-slm.R`: Analysis of the Boston housing data set using the main spatial econometrics models.
* `boston-slm-impacts.R`: Computation of the impacts for the Boston housing data example.
* `boston-slm-full.R`: Analysis of the Boston housing data set using the main spatial econo- metrics models and the full adjacency matrix to perform prediction on the missing values.
* `katrina-slm.R`: Analysis of the Katrina business data using the main spatial econometrics models with a spatial probit.
* `katrina-slm-neigh.R`: Selection of the number of optimal nearest neighbours for the adjacency matrix using the Katrina business data.
* `katrina-slm-impacts.R`: Computation of the impacts for the Katrina business data example.
* `katrina-spatialprobit.R`: Code to fit the spatial probit models (and compute their impacts) using MCMC with the `spatialprobit` package.

In addition, we provide a customized version of package `spatialprobit` (`spatialprobit_0.9-11.tar.gz`) that sets the variacne of the random effect to 1 in the SEM probit model. It is already set to 1 in the SAR probit model.



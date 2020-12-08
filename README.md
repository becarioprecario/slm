# Spatial lag model with INLA: the `slm` latent effect

This repository contains the `R` code of paper [*Estimating Spatial Econometrics Models with Integrated Nested Laplace Approximation*](https://arxiv.org/abs/1703.01273) by Virgilio Gómez-Rubio, Roger S. Bivand, Håvard Rue.


Files available are:

* `boston-slm.R`: Analysis of the Boston housing data set using the main spatial econometrics models.
* `boston-slm-impacts.R`: Computation of the impacts for the Boston housing data example.
* `boston-slm-full.R`: Analysis of the Boston housing data set using the main spatial econo- metrics models and the full adjacency matrix to perform prediction on the missing values.
* `katrina-slm.R`: Analysis of the Katrina business data using the main spatial econometrics models with a spatial probit.
* `katrina-slm-neigh.R`: Selection of the number of optimal nearest neighbours for the adjacency matrix using the Katrina business data.
* `katrina-slm-impacts.R`: Computation of the impacts for the Katrina business data example.


**Note that this is old code and that may not work straight away with the latest version of some spatial packages.**

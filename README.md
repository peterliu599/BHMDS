# Bayesian Hyperbolic Multidimensional Scaling (BHMDS)
By Bolun (Peter) Liu, Department of Biostatistics, Johns Hopkins Bloomberg School of Public Health
Shane Lubold, Adrian E. Raftery, Tyler H. McCormick, Department of Statistics, University of Washington

This repository contains the code to reproduce the figures in the paper "Bayesian Hyperbolic Multidimensional Scaling".

URL: https://arxiv.org/abs/2210.15081. 

There are four folders in this repository:
1) "simulation": Contains code for reproducing Figure 1-4, i.e. simulation experiments and calibration evaluation for BHMDS;
2) "stress & distortion": Contains code for reproducing the stress and distortion evaluation in Table 1 & 2 for data sets
- Karate club
- Phylogenetic tree
- CS phd
- Wordnet mammal subtree
We provide the original data sets as well as bhmds & bmds embeddings.
3) "log likelihood comparison wordnet": Contains code for reproducing Figure 5, i.e. the log-likelihood comparison between full & case-control approximated MCMC using Wordnet mammal subtree data sets. 
4) "application": Contain code for reproducing Figure 6 & 7, i.e. the cluster-wise cell type distance and frequency of their rank statistics. 

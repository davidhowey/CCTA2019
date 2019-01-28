# CCTA2019 Readme

Matlab code for Paper submitted to 
the Third IEEE Conference on Control Technology and Applications (CCTA 2019).

Required: Mosek-SDP package and Yalmip
crossover_data.mat: experimental data in Matlab format
crossover_lmi_polytopic_design.m: computation of the observer gain via polytopic solution of the LMI problem--> obs_gains_poly.mat
crossover_sys_data.m: simulation of battery system with linear crossover dynamic and comparison with experimental data
crossover_obs_data.m: Simulation of Augmented State Observer, battery system with linear crossover dynamic and experimetal data -->crossover_obs_data_poly.mat

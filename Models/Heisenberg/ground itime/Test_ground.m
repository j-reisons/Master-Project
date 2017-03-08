clc
close all
clear all
%%
J = 1;
U = 1;
N = 10;
D = 20;

H_MPO = Heisenberg_H(N,J,U);

ground_MPS = Heisenberg_ground_itime(N,U,J,D);


%%
(braket(ground_MPS,apply(H_MPO,apply(H_MPO,ground_MPS))) - braket(ground_MPS,apply(H_MPO,ground_MPS))^2)


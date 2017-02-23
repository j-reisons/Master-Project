clc
close all
clear all
%%
D_init = 50;
D_comp = 10;
d = 2;
N = 50;
State = random_mps(N,D_init,d);
%% SVD compression as initial guess
State_SVD = State;
State_SVD = sweep(State_SVD,1);
State_SVD = sweep(State_SVD,-1,D_comp);
State_Iter = Iter_comp(State,D_comp,1E-17);

%%
clc
close all
clear all

addpath(genpath('C:\Users\Juris\Documents\MATLAB\Master Project\General functions'))
%%
fprintf('Generating random MPS... ');
D_init = 250;
D_comp = 50;
d = 2;
N = 50;
tol = 1E-12;
State = random_mps(N,D_init,d);
fprintf('done \n');
%% 
State_SVD = State;
fprintf('Testing SVD compression... ');
tic
State_SVD = sweep(State_SVD,1);
State_SVD = sweep(State_SVD,-1,tol,D_comp);

fprintf('done. \nAccuracy is  ')
fprintf(num2str(norm(1 - braket(State,State_SVD))))
fprintf('\n')
toc

fprintf('Testing iterative compression... ');
tic
State_Iter = Iter_comp(State,D_comp,tol);

fprintf('done. \nAccuracy is ')
fprintf(num2str(norm(1 - braket(State,State_Iter))))
fprintf('\n')
toc

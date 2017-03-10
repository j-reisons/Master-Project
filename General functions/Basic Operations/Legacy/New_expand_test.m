clc
close all
clear all
%%
A = random_mps(12,100,2);

fprintf('Old took\n')
tic
old = expand_MPS_legacy(A);
toc


fprintf('New took\n')
tic
new = expand_MPS(A);
toc
%%
assert(approx(norm(old-new),0,1E-10));
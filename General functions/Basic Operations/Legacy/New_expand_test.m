clc
close all
clear all
%%
A = random_mps(12,100,2);

old = expand_MPS_legacy(A);
new = expand_MPS(A);

%%
assert(approx(norm(old-new),0,1E-10));
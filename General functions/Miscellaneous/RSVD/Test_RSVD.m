%RSVD test
clc
close all
clear all
%%
A = rand(1600,400) - 0.5;

tic
[U,S,V] = svd(A,'econ');
toc

tic
[U2,S2,V2] = rsvd(A,100,5);
toc

%%

S = diag(S);
S2 = diag(S2);

plot(S)
hold on
plot(S2);
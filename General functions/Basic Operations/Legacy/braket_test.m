clc
close all
clear all

%%
N = 20;
D_1 = 50;
D_2 = 15;
d = 2;

%%
A = random_mps(N,D_1,d);
B = random_mps(N,D_2,d);


fprintf('Current version took \n');
tic
braket1 = braket(A,B);
toc

fprintf('worse version took \n');
tic
braket2 = worse_braket(A,B);
toc

fprintf('even worse version took \n');
tic
braket3 = even_worse_braket(A,B);
toc

%%

assert(approx(braket1,braket2,1E-10));
assert(approx(braket1,braket3,1E-10));
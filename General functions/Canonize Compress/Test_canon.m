clc
clear all
close all
%%
N = 15;
D = 20;
d = 2;

fprintf('Testing canonization functions...')
A = random_mps(15,20,2);
assert(isequal(iscanon(A),1));
%%
A = sweep(A,-1);
assert(isequal(iscanon(A),-1));
%%
O = Heisenberg_H(N,1,1);
A = apply(O,A);
assert(isequal(iscanon(A),0));

fprintf('done.\n')
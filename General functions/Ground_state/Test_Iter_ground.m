clc
close all
clear all
%%

N = 40;
J = 1;
U = 1;
D = 50;
tol = 1E-9;

%%
Ham = Heisenberg_H(N,J,U);

fprintf('Generating random MPS...')
MPS = random_mps(N,D,2,-1);
fprintf('done.\n')

fprintf('Testing ground state search...')

[MPS,energy] = Iter_ground(Ham,MPS,tol);

%%

HMPS = apply(Ham,MPS);

HH = braket(HMPS,HMPS);
H_H = braket(MPS,HMPS)^2;
test = (HH - H_H)/H_H;

assert(approx(test,0,tol));
fprintf('done.\n')
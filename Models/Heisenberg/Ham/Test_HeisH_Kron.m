clc
close all
clear all

%%
fprintf('Testing MPO form of XXZ Hamiltonian...')
N = 8;
J = 1;
U = 2;
d = 2;

%%
MPO = Heisenberg_H(N,J,U);

S_X =[0,1;1,0];
S_Y = [0,-1i;1i,0];
S_Z = [1,0;0,-1];
Ham_pair = J*(kron(S_X,S_X) + kron(S_Y,S_Y)) + U*(kron(S_Z,S_Z));

H_kron = 0;

for i = 1:1:N-1
    H_kron = H_kron + kron(kron(eye(d^(i-1)),Ham_pair),eye(d^(N-i-1)));
end

H_MPO = expand_MPO(MPO);


assert(isequal(H_kron,H_MPO));

fprintf(' done.\n')
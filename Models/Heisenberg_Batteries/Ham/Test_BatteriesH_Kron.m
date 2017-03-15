clc
close all
clear all
%%
fprintf('Testing MPO form of Battery Hamiltonian...')
N = 4;
J = 1;
Uc = 1;
Ub = 4;
d = 2;
dim = 2^(3*N);

%%
MPO = Batteries_H(N,J,Uc,Ub);
%% Kroning

S_X =[0,1;1,0];
S_Y = [0,-1i;1i,0];
S_Z = [1,0;0,-1];

Ham_pair_c = sparse(J*(kron(S_X,S_X) + kron(S_Y,S_Y)) + Uc*(kron(S_Z,S_Z)));
Ham_pair_b = sparse(J*(kron(S_X,S_X) + kron(S_Y,S_Y)) + Ub*(kron(S_Z,S_Z)));

H_kron = sparse(dim,dim);
for i = 1:1:N
    H_kron = H_kron + kron(kron(speye(d^(i-1)),Ham_pair_b),speye(d^((3*N)-i-1)));
end
for i = N+1:1:(2*N)-1
    H_kron = H_kron + kron(kron(speye(d^(i-1)),Ham_pair_c),speye(d^((3*N)-i-1)));
end
for i = 2*N:1:(3*N)-1
    H_kron = H_kron + kron(kron(speye(d^(i-1)),Ham_pair_b),speye(d^((3*N)-i-1)));
end

%%
H_MPO = expand_MPO(MPO);

assert(approx(max(max(H_kron-H_MPO)),0,1E-15));
max(max(H_kron-H_MPO));

fprintf(' done.\n')
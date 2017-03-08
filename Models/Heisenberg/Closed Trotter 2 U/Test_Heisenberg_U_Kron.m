close all
clear all
clc

addpath('General functions','Models/Heisenberg_Batteries')

%%

N = 9;
J = 2;
U = 1;
d = 2;
dt = 0.01; % Time increment

[U_even_dt,U_odd_dt] = Heisenberg_U(N,J,U,dt);
[U_even_half,U_odd_half] = Heisenberg_U(N,J,U,dt/2);
U_total = compressMPO(U_odd_half,U_even_dt,U_odd_half);

%% Kroning
fprintf('Generating Kroneker products...');

S_X =[0,1;1,0];
S_Y = [0,-1i;1i,0];
S_Z = [1,0;0,-1];
Ham_pair = J*(kron(S_X,S_X) + kron(S_Y,S_Y)) + U*(kron(S_Z,S_Z));

H_kron = 0;
H_even_kron = 0;
H_odd_kron = 0;

for i = 1:1:N-1
    H_kron = H_kron + kron(kron(eye(d^(i-1)),Ham_pair),eye(d^(N-i-1)));
end

for i = 1:2:N-1
    H_odd_kron = H_odd_kron + kron(kron(eye(d^(i-1)),Ham_pair),eye(d^(N-i-1)));
end

for i = 2:2:N-1
    H_even_kron = H_even_kron + kron(kron(eye(d^(i-1)),Ham_pair),eye(d^(N-i-1)));
end

U_even_kron = expm(-1i*dt*H_even_kron);
U_odd_kron = expm(-1i*dt*H_odd_kron);
U_kron = expm(-1i*dt*H_kron);
Trotter_kron = expm(-1i*(dt/2)*H_odd_kron)*expm(-1i*dt*H_even_kron)*expm(-1i*(dt/2)*H_odd_kron);

norm(Trotter_kron - U_kron,'fro')

fprintf('\tdone\n');
%% MPS U
fprintf('Expanding MPOs...');

Test_even = expand_MPO(U_even_dt);
Test_odd = expand_MPO(U_odd_dt);
Test = expand_MPO(U_total);

fprintf('\tdone\n');
%%
tolerance = 1e-10;

fprintf('Testing U_even...');
assert(approx(norm(U_even_kron - Test_even,'fro'),0,tolerance));
fprintf('\tdone\n');

fprintf('Testing U_odd...');
assert(approx(norm(U_odd_kron - Test_odd,'fro'),0,tolerance));
fprintf('\tdone\n');

fprintf('Testing U_trotter...');
assert(approx(norm(Trotter_kron - Test,'fro'),0,tolerance));
fprintf('\tdone\n');

fprintf('Testing U_norm...');
assert(approx(norm(Test,'fro'),2^(N/2),tolerance));
fprintf('\tdone\n');

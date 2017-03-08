close all
clear all
clc

addpath('General functions','Models/Heisenberg_Batteries')

%%

N = 2;
noteven = mod(N,2);
J = 2;
U_b = 0;
U_c = 1;
d = 2;
dt = 0.02; % Time increment

[U_even_dt,U_odd_dt] = Heisenberg_Batteries_U(N,J,U_b,U_c,dt);
[U_even_half,U_odd_half] = Heisenberg_Batteries_U(N,J,U_b,U_c,dt/2);
U_total = compressMPO(U_odd_half,U_even_dt,U_odd_half);

%% Kroning
fprintf('Generating Kroneker products...');

S_X =[0,1;1,0];
S_Y = [0,-1i;1i,0];
S_Z = [1,0;0,-1];


Ham_pair_c = J*(kron(S_X,S_X) + kron(S_Y,S_Y)) + U_c*(kron(S_Z,S_Z));
Ham_pair_b = J*(kron(S_X,S_X) + kron(S_Y,S_Y)) + U_b*(kron(S_Z,S_Z));

H_kron = 0;
H_even_kron = 0;
H_odd_kron = 0;


for i = 1:1:N
    H_kron = H_kron + kron(kron(eye(d^(i-1)),Ham_pair_b),eye(d^((3*N)-i-1)));
end
for i = N+1:1:(2*N)-1
    H_kron = H_kron + kron(kron(eye(d^(i-1)),Ham_pair_c),eye(d^((3*N)-i-1)));
end
for i = 2*N:1:(3*N)-1
    H_kron = H_kron + kron(kron(eye(d^(i-1)),Ham_pair_b),eye(d^((3*N)-i-1)));
end

for i = 1:2:N
    H_odd_kron = H_odd_kron + kron(kron(eye(d^(i-1)),Ham_pair_b),eye(d^((3*N)-i-1)));
end
for i = 2:2:N
    H_even_kron = H_even_kron + kron(kron(eye(d^(i-1)),Ham_pair_b),eye(d^((3*N)-i-1)));
end

if noteven
    for i = N+1:2:(2*N)-1
        H_even_kron = H_even_kron + kron(kron(eye(d^(i-1)),Ham_pair_c),eye(d^((3*N)-i-1)));
    end
    
    for i = N+2:2:(2*N)-1
        H_odd_kron = H_odd_kron + kron(kron(eye(d^(i-1)),Ham_pair_c),eye(d^((3*N)-i-1)));
    end
else
    for i = N+1:2:(2*N)-1
        H_odd_kron = H_odd_kron + kron(kron(eye(d^(i-1)),Ham_pair_c),eye(d^((3*N)-i-1)));
    end  
    for i = N+2:2:(2*N)-1
        H_even_kron = H_even_kron + kron(kron(eye(d^(i-1)),Ham_pair_c),eye(d^((3*N)-i-1)));
    end
end

for i = 2*N:2:3*N-1
    H_even_kron = H_even_kron + kron(kron(eye(d^(i-1)),Ham_pair_b),eye(d^((3*N)-i-1)));
end
for i = 2*N+1:2:3*N-1
    H_odd_kron = H_odd_kron + kron(kron(eye(d^(i-1)),Ham_pair_b),eye(d^((3*N)-i-1)));
end

U_even_kron = expm(-1i*dt*H_even_kron);
U_odd_kron = expm(-1i*dt*H_odd_kron);
U_kron = expm(-1i*dt*H_kron);
Trotter_kron = expm(-1i*(dt/2)*H_odd_kron)*expm(-1i*dt*H_even_kron)*expm(-1i*(dt/2)*H_odd_kron);

fprintf('\tdone\n');
%% MPS U
fprintf('Expanding MPOs...');

Test_even = expand_MPO(U_even_dt);
Test_odd = expand_MPO(U_odd_dt);
Test = expand_MPO(U_total);

fprintf('\tdone\n');
%%
nazi_tolerance = 1e-12;
compressed_tolerance = 1e-5;

fprintf('Testing U_even...');
assert(approx(norm(U_even_kron - Test_even,'fro'),0,nazi_tolerance));
fprintf('\tdone\n');

fprintf('Testing U_odd...');
assert(approx(norm(U_odd_kron - Test_odd,'fro'),0,nazi_tolerance));
fprintf('\tdone\n');

fprintf('Testing U_trotter...');
%norm(Trotter_kron - Test,'fro')
assert(approx(norm(Trotter_kron - Test,'fro'),0,compressed_tolerance));
fprintf('\tdone\n');

fprintf('Testing U_norm...');
assert(approx(norm(Test,'fro'),2^(N*(3/2)),nazi_tolerance));
fprintf('\tdone\n');

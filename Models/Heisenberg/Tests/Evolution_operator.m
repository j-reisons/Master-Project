close all
clear all
clc
%%

N = 8;
J = 1;
U = 1;
d = 2;
dt = 0.1; % Time increment

[U_odd,U_even] = Heisenberg_U(N,J,U,dt);
%% Stuff
S_X =[0,1;1,0];
S_Y = [0,-1i;1i,0];
S_Z = [1,0;0,-1];
Ham_pair = J*(kron(S_X,S_X) + kron(S_Y,S_Y)) + U*(kron(S_Z,S_Z));
U_pair = expm(-1i*dt*Ham_pair);
%% Kron
 H_odd_kron = 0;

 for i = 1:2:N-1
    H_odd_kron = H_odd_kron + kron(kron(eye(d^(i-1)),Ham_pair),eye(d^(N-i-1)));
 end
 
 U_odd_kron = expm(-1i*dt*H_odd_kron);
 
H_even_kron = 0;
for i = 2:2:N-1
    H_even_kron = H_even_kron + kron(kron(eye(d^(i-1)),Ham_pair),eye(d^(N-i-1)));
end

U_even_kron = expm(-1i*dt*H_even_kron);

%% MPS U
Test_odd = expand_MPO(U_odd);
Test_even = expand_MPO(U_even);
%%
norm(Test_odd - U_odd_kron)
norm(Test_even - U_even_kron)
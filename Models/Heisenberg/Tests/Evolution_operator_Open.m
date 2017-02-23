close all
clear all
clc
%%
N = 5;
J = 1;
U = 1;
G = 1;
d = 2;
dt = 0.05; % Time increment

[U_odd,U_even] = HeisenbergOpen_U_O2(N,J,U,G,dt); %HeisenbergOpenDisordered(N,J,U,G,0,dt);
%% Hamiltonian
S_X =[0,1;1,0];
S_Y = [0,-1i;1i,0];
S_Z = [1,0;0,-1];
S_plus = [0,1;0,0];
S_minus = S_plus.';
I = eye(2);
S_plus_1 = kron(S_plus,speye(2^(N-1)));
S_minus_N = kron(speye(2^(N-1)),S_minus);

H = zeros(2^N);
for i = 1:N-1
    O_X = 1;
    O_Y = 1;
    O_Z = 1;
    for j = 1:N
        if j == i || j == i+1
            O_X = kron(O_X,S_X);
            O_Y = kron(O_Y,S_Y);
            O_Z = kron(O_Z,S_Z);
        else
            O_X = kron(O_X,I);
            O_Y = kron(O_Y,I);
            O_Z = kron(O_Z,I);
        end
    end
    H = H + J*O_X + J*O_Y + U*O_Z ;
end

H_eff = H -(((1i)*G)/2)*(S_plus_1'*S_plus_1) - (((1i)*G)/2)*(S_minus_N'*S_minus_N);
U_eff = expm(-1i*H_eff*dt);

Ham_pair = J*(kron(S_X,S_X) + kron(S_Y,S_Y)) + U*(kron(S_Z,S_Z));
H_odd_kron = 0;
 for i = 1:2:N-1
    H_odd_kron = H_odd_kron + kron(kron(eye(d^(i-1)),Ham_pair),eye(d^(N-i-1)));
 end
H_odd_kron = H_odd_kron -(((1i)*G)/4)*(S_plus_1'*S_plus_1) - (((1i)*G)/4)*(S_minus_N'*S_minus_N);
U_odd_kron = expm(-1i*(dt/2)*H_odd_kron);
H_even_kron = 0;
for i = 2:2:N-1
    H_even_kron = H_even_kron + kron(kron(eye(d^(i-1)),Ham_pair),eye(d^(N-i-1)));
end
H_even_kron = H_even_kron -(((1i)*G)/4)*(S_plus_1'*S_plus_1) - (((1i)*G)/4)*(S_minus_N'*S_minus_N);
U_even_kron = expm(-1i*dt*H_even_kron);

%%
Test_odd = expand_MPO(U_odd);
Test_even = expand_MPO(U_even);
Test = Test_odd*Test_even*Test_odd;
%%
norm(Test_odd - U_odd_kron)
norm(Test_even - U_even_kron)
norm(Test - U_eff)

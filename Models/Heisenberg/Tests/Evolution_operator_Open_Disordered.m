close all
clear all
clc
%%
N = 8;
J = 1;
U = 1;
G = 1;
D = 0;
d = 2;
dt = 0.05; % Time increment

[U_L,U_odd,U_even,disorder] = HeisenbergOpenDisordered_H_L_trotter(N,J,U,G,D,dt);
[U_odd_2,U_even_2] = HeisenbergOpenDisordered(disorder,N,J,U,G,D,dt);
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
for i=1:N
    H_eff = H_eff + D*disorder(i)*kron(kron(eye(2^(i-1)),S_Z),eye(2^(N-i)));
end
U_kron = expm(-1i*dt*H_eff);
%%
Test_L = expand_MPO(U_L);
Test_odd = expand_MPO(U_odd);
Test_even = expand_MPO(U_even);
%%
%norm(Test_odd - U_odd_kron)
%norm(Test_even - U_even_kron)
norm(Test_L*Test_odd*Test_even*Test_odd*Test_L - U_kron)
%%
Test_odd_2 = expand_MPO(U_odd_2);
Test_even_2 = expand_MPO(U_even_2);
%%
%norm(Test_odd - U_odd_kron)
%norm(Test_even - U_even_kron)
norm(Test_odd_2*Test_even_2*Test_odd_2 - U_kron)
%%
sum(sum(abs(Test_L*Test_odd*Test_even*Test_odd*Test_L - Test_odd_2*Test_even_2*Test_odd_2)))
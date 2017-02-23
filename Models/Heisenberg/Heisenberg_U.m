function [U_odd,U_even] = Heisenberg_U(N,J,U,dt)
%HEISENBERG Returns odd and even evolution operators for the Heisenberg
%Hamiltonian of parameters N,J,U,dt
%Outlined in Schollwock 7.1.1 p75 and 7.1.2 p 77
%% Pauli and co
d = 2;

S_X =[0,1;1,0];

S_Y = [0,-1i;1i,0];

S_Z = [1,0;0,-1];
%% Pair evolution operator

Ham_pair = J*(kron(S_X,S_X) + kron(S_Y,S_Y)) + U*(kron(S_Z,S_Z));
U_pair = expm(-1i*dt*Ham_pair);

%% U -> P reshaping
P = reshape(U_pair,[2,2,2,2]);
P = permute(P,[4,2,3,1]);
P = reshape(P,[d^2,d^2]); % (sig 1 sig 1'),(sig2 sig2')

%% SVD of P

[U,S,V] = svd(P,'econ');
k = size(S,2);
U = U*sqrt(S); %(sig1 sig1'), k
U_bar = sqrt(S)*V'; %k,(sig2 sig2')


%% U, U_bar and eye(2) reshaping

U = reshape(U,[d,d,1,k]);% sig 1,sig 1',1,k
U = permute(U,[3,4,1,2]);%1,k,sig1,sig1'
U_bar = reshape(U_bar,[k,1,d,d]);% k, 1 ,sig2 sig2'
I_site = reshape(eye(2),[1,1,d,d]);% 1,1,sig,sig'
%% Putting it all into U_even and U_odd

U_odd = cell(1,N);

for i = 1:2:N-1;
U_odd{i}= U;
U_odd{i+1}= U_bar;
end

U_even = cell(1,N);

U_even{1}= I_site;

for i = 2:2:N-1;
U_even{i} = U;
U_even{i+1} = U_bar;
end

% Chain length problems
if not(mod(N,2))
    U_even{N} = I_site;
else
    U_odd{N} = I_site;
end

end
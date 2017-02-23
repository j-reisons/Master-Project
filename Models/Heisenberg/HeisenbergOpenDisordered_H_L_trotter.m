function [U_L,U_odd,U_even,disorder] = HeisenbergOpenDisordered_H_L_trotter(N,J,U,G,D,dt)
%HEISENBERGOPEN_U Returns second order Trotter decomposition
%of non unitary evolution operators for the
%Heisenberg Hamiltonian of parameters N,J,U,G,dt
% U_odd(dt/2) and U_even(dt); U(dt) ~= U_odd(dt/2)U_even(dt)U_odd(dt/2)

%% Pauli and co
d = 2;

S_X =[0,1;1,0];

S_Y = [0,-1i;1i,0];

S_Z = [1,0;0,-1];

S_plus = [0,1;0,0];

S_minus = S_plus.';

disorder = rand(N,1) - 0.5;
%% Pair evolution operators

U_cell = cell(2,1);
U_bar_cell = cell(2,1);

Ham_pair = J*(kron(S_X,S_X) + kron(S_Y,S_Y)) + U*(kron(S_Z,S_Z));

for i = 1:2
    %% U -> P reshaping
    U_pair = expm(-1i*(dt/i)*Ham_pair);
    P = reshape(U_pair,[2,2,2,2]);
    P = permute(P,[4,2,3,1]); %
    P = reshape(P,[d^2,d^2]); % (sig 1 sig 1'),(sig2 sig2')
    
    %% SVD of P
    [U,S,V] = svd(P,'econ');
    k = size(S,2);
    U = U*sqrt(S); %(sig1 sig1'), k
    U_bar = sqrt(S)*V'; %k,(sig2 sig2')

    %% U, U_bar and eye(2) reshaping    
    U = reshape(U,[d,d,1,k]);% sig 1,sig 1',1,k
    U_cell{i} = permute(U,[3,4,1,2]);%1,k,sig1,sig1'
    U_bar_cell{i} = reshape(U_bar,[k,1,d,d]);% k, 1 ,sig2 sig2'
    
end
I_site = reshape(eye(2),[1,1,d,d]);% 1,1,sig,sig'

U_odd = cell(1,N);

for i = 1:2:N-1;
    U_odd{i}= U_cell{2};
    U_odd{i+1}= U_bar_cell{2};
end

U_even = cell(1,N);
U_even{1}= I_site;
for i = 2:2:N-1;
    U_even{i} = U_cell{1};
    U_even{i+1} = U_bar_cell{1};
end

% Chain length problems
if not(mod(N,2))
    U_even{N} = I_site;
else
    U_odd{N} = I_site;
end

%% U_L
L = cell(1,N);

L{1} = (-1i)*G*0.5*S_minus*S_plus + D*disorder(1)*S_Z;
L{N} = (-1i)*G*0.5*S_plus*S_minus + D*disorder(N)*S_Z;
for i = 2:N-1
    L{i} = D*disorder(i)*S_Z;
end

U_L = cell(1,N);
for i = 1:N
U_L{i} = expm(-1i*dt*0.5*L{i});
U_L{i} = reshape(U_L{i},[1,1,d,d]);
end

end
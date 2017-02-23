function [U_odd_1_2,U_even_1_2,U_odd_3,U_even_3] = HeisenbergOpen_U_O4(N,J,U,G,dt)
%HEISENBERGOPEN_U Returns second order Trotter decomposition 
%of non unitary evolution operators for the 
%Heisenberg Hamiltonian of parameters N,J,U,G,dt
% U_odd(dt/2) and U_even(dt); U(dt) ~= U_odd(dt/2)U_even(dt)U_odd(dt/2)
% 
%% Pauli and co
d = 2;

S_X =[0,1;1,0];

S_Y = [0,-1i;1i,0];

S_Z = [1,0;0,-1];

S_plus = [0,1;0,0];

S_minus = S_plus.';

tau1_2 = dt/(4 - 4^(1/3));
tau3 = dt - 4*tau1_2; 

tau = [tau1_2,tau1_2/2,tau3,tau3/2];
%% Pair evolution operator, for middle, first and last bonds
Ham_pair = {};

Ham_pair{1} = J*(kron(S_X,S_X) + kron(S_Y,S_Y)) + U*(kron(S_Z,S_Z));
Ham_pair{2} = Ham_pair{1} - (1i)*(G/4)*kron(S_plus'*S_plus,eye(2));
Ham_pair{3} = Ham_pair{1} - (1i)*(G/4)*kron(eye(2),S_minus'*S_minus);

U_cell = cell(3,2);
U_bar_cell = cell(3,2);
for l = 1:3 % Middle, first and last site U matrices
    for j = 1:length(tau) % various tau lengths
        U_pair = expm(-1i*(tau(j)*Ham_pair{l})); % 2 site evolution operator
        %% U -> P reshaping
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
        U = permute(U,[3,4,1,2]);%1,k,sig1,sig1'
        U_bar = reshape(U_bar,[k,1,d,d]);% k, 1 ,sig2 sig2'
        
        U_cell{l,j} = U;
        U_bar_cell{l,j} = U_bar;
    end
end
%% Putting it all into U_even and U_odd

U_odd_1_2 = cell(1,N);
U_odd_1_2{1} = U_cell{2,2};
U_odd_1_2{2} = U_bar_cell{2,2};

U_odd_3 = cell(1,N);
U_odd_3{1} = U_cell{2,4};
U_odd_3{2} = U_bar_cell{2,4};

for i = 3:2:N-1;
U_odd_1_2{i}= U_cell{1,2};
U_odd_1_2{i+1}= U_bar_cell{1,2};
U_odd_3{i} = U_cell{1,4};
U_odd_3{i+1} = U_bar_cell{1,4};
end

U_even_1_2 = cell(1,N);
U_even_1_2{1}= reshape(expm(-(G/4)*tau(1)*(S_plus')*S_plus),[1,1,d,d]);

U_even_3 = cell(1,N);
U_even_3{1}= reshape(expm(-(G/4)*tau(3)*(S_plus')*S_plus),[1,1,d,d]);

for i = 2:2:N-1;
U_even_1_2{i} = U_cell{1,1};
U_even_1_2{i+1} = U_bar_cell{1,1};
U_even_3{i} = U_cell{1,3};
U_even_3{i+1} = U_bar_cell{1,3};
end

% Chain length matters
if not(mod(N,2)) %even
    U_odd_1_2{N-1}= U_cell{3,2};
    U_odd_1_2{N} = U_bar_cell{3,2};
    U_odd_3{N-1}= U_cell{3,4};
    U_odd_3{N} = U_bar_cell{3,4};
    
    U_even_1_2{N} = reshape(expm(-(G/4)*tau(1)*(S_minus')*S_minus),[1,1,d,d]);
    U_even_3{N} = reshape(expm(-(G/4)*tau(3)*(S_minus')*S_minus),[1,1,d,d]);
else %odd
    U_even_1_2{N-1} = U_cell{3,1};
    U_even_1_2{N} = U_bar_cell{3,1};
    U_even_3{N-1} = U_cell{3,3};
    U_even_3{N} = U_bar_cell{3,3};
    
    U_odd_1_2{N} = reshape(expm(-(G/4)*tau(2)*(S_minus')*S_minus),[1,1,d,d]);
    U_odd_3{N} = reshape(expm(-(G/4)*tau(4)*(S_minus')*S_minus),[1,1,d,d]);
end

end
function [U_odd,U_even] = HeisenbergOpenDisordered(N,J,U,G,D,dt)
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
%% Pair evolution operators for all bonds
Base_ham_pair = J*(kron(S_X,S_X) + kron(S_Y,S_Y)) + U*(kron(S_Z,S_Z));
Ham_pair = cell(N-1,1);

Ham_pair{1} = Base_ham_pair - (1i)*(G/4)*kron(S_plus'*S_plus,eye(2)) ...
    + 0.5*D*(disorder(1)*kron(S_Z,eye(2))+ disorder(2)*kron(eye(2),S_Z));

for i = 2:N-2
    Ham_pair{i} = Base_ham_pair + 0.5*D*(disorder(i)*kron(S_Z,eye(2)) ...
        + disorder(i+1)*kron(eye(2),S_Z));
end

Ham_pair{N-1} = Base_ham_pair -(1i)*(G/4)*kron(eye(2),S_minus'*S_minus) ...
 + 0.5*D*(disorder(N-1)*kron(S_Z,eye(2))+ disorder(N)*kron(eye(2),S_Z));

U_cell = cell(N-1,2);
U_bar_cell = cell(N-1,2);

for l = 1:N-1 % Over bonds
    for j = [1,2] % dt and dt/2 for even and odd operators
        U_pair = expm(-1i*(dt/j)*Ham_pair{l}); % 2 site evolution operator
        %% U -> P reshaping
        P = reshape(U_pair,[2,2,2,2]);
        % P = permute(P,[1,3,2,4]);
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
        
        U_cell{l,j} = U;
        U_bar_cell{l,j} = U_bar;
    end
end
%% Putting it all into U_even and U_odd

U_odd = cell(1,N);

for i = 1:2:N-1;
U_odd{i}= U_cell{i,2};
U_odd{i+1}= U_bar_cell{i,2};
end

U_even = cell(1,N);
U_even{1}= reshape(expm((-1i*dt)*(-1i*(G/4)*(S_plus')*S_plus + ...
    0.5*D*disorder(1)*S_Z)),[1,1,d,d]);

for i = 2:2:N-1;
U_even{i} = U_cell{i,1};
U_even{i+1} = U_bar_cell{i,1};
end

% Chain length matters
if not(mod(N,2)) %even
    U_odd{N-1}= U_cell{N-1,2};
    U_odd{N} = U_bar_cell{N-1,2};
    U_even{N} = reshape(expm((-1i*dt)*(-1i*(G/4)*(S_minus')*S_minus + ...
        0.5*D*disorder(N)*S_Z)),[1,1,d,d]);
else %odd
    U_even{N-1} = U_cell{N-1,1};
    U_even{N} = U_bar_cell{N-1,1};
    U_odd{N} = reshape(expm((-1i*dt/2)*(-1i*(G/4)*(S_minus')*S_minus + ...
        0.5*D*disorder(N)*S_Z)),[1,1,d,d]);
end

end
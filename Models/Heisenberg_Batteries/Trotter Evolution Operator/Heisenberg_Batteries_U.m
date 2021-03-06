function [U_even,U_odd] = Heisenberg_Batteries_U(Nc,Nb,J,U_b,U_c,dt)
%HEISENBERG_Batteries Returns Even and odd evolution operators
%Hamiltonian of parameters N,J,U,dt
%Outlined in Schollwock 7.1.1 p75 and 7.1.2 p 77
%% Pauli and co
d = 2;

notevenB = mod(Nb,2);
notevenBC = mod(Nb+Nc,2);
notevenBCB = mod(2*Nb+Nc,2);

S_X =[0,1;1,0];

S_Y = [0,-1i;1i,0];

S_Z = [1,0;0,-1];

U_parameters = cell(1,2);
U_parameters{1} = U_b;
U_parameters{2} = U_c;

U_c = cell(1,2);
U_bar_c = cell(1,2);

%% Pair evolution operator

for i = 1:2
    Ham_pair = J*(kron(S_X,S_X) + kron(S_Y,S_Y)) + U_parameters{i}*(kron(S_Z,S_Z));
    U_pair = expm(-1i*dt*Ham_pair);
    
    % U -> P reshaping
    P = reshape(U_pair,[2,2,2,2]);
    P = permute(P,[4,2,3,1]);
    P = reshape(P,[d^2,d^2]); % (sig 1 sig 1'),(sig2 sig2')
    
    % SVD of P
    
    [U,S,V] = svd(P,'econ');
    k = size(S,2);
    U = U*sqrt(S); %(sig1 sig1'), k
    U_bar = sqrt(S)*V'; %k,(sig2 sig2')
    
    
    % U, U_bar and eye(2) reshaping
    
    U = reshape(U,[d,d,1,k]);% sig 1,sig 1',1,k
    U_c{i} = permute(U,[3,4,1,2]);%1,k,sig1,sig1'
    U_bar_c{i} = reshape(U_bar,[k,1,d,d]);% k, 1 ,sig2 sig2'
end

I_site = reshape(eye(2),[1,1,d,d]);% 1,1,sig,sig'

%% Odd bonds

U_odd = cell(1,Nc + 2*Nb);

for i = 1 : 2 : Nb;
    U_odd{i}= U_c{1};
    U_odd{i+1}= U_bar_c{1};
end

if notevenB
    for i = Nb+2 : 2 : Nb + Nc - 1;
    U_odd{i}= U_c{2};
    U_odd{i+1}= U_bar_c{2};
    end
else
    for i = Nb+1 : 2 : Nb + Nc - 1;
    U_odd{i}= U_c{2};
    U_odd{i+1}= U_bar_c{2};
    end
end

if notevenBC
    for i = Nb + Nc : 2 : 2*Nb + Nc - 1;
        U_odd{i}= U_c{1};
        U_odd{i+1}= U_bar_c{1};
    end
else
    for i = Nb + Nc + 1 : 2 : 2*Nb + Nc - 1;
        U_odd{i}= U_c{1};
        U_odd{i+1}= U_bar_c{1};
    end
end

%Chain length parity
if notevenBCB
    U_odd{Nb*2 + Nc} = I_site;
end


%% Even bonds

U_even = cell(1,Nc + 2*Nb);
U_even{1}= I_site;

for i = 2 : 2 : Nb;
    U_even{i} = U_c{1};
    U_even{i+1} = U_bar_c{1};
end

if notevenB
    for i = Nb+1 : 2 :Nb + Nc -1;
        U_even{i} = U_c{2};
        U_even{i+1} = U_bar_c{2};
    end   
else
    for i = Nb+2 : 2 :Nb + Nc -1;
        U_even{i} = U_c{2};
        U_even{i+1} = U_bar_c{2};
    end
end

if notevenBC
    for i = Nb + Nc + 1 : 2 : 2*Nb + Nc -1
        U_even{i} = U_c{1};
        U_even{i+1} = U_bar_c{1};
    end
else
    for i = Nb + Nc : 2 : 2*Nb + Nc -1
        U_even{i} = U_c{1};
        U_even{i+1} = U_bar_c{1};
    end
end

%Chain length parity
if ~notevenBCB
    U_even{Nb*2 + Nc} = I_site;
end

end


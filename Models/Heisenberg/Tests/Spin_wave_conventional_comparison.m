clc
close all
clear all
%% Parameters
N = 9;
J = 1;
U = 0;
T = 10;
dt = 5*10^-2;
time = 0:dt:T;
%% Pauli and co.
S_X = sparse([
    [0 , 1]
    [1 , 0]
    ]);

S_Y = sparse([
    [0 , -1i]
    [1i , 0]
    ]);

S_Z = sparse([
    [1 , 0]
    [0 ,-1]
    ]);

S_plus = sparse([
    [0,1]
    [0,0]
    ]);

S_minus = S_plus.';

I = speye(2);

%% Building the Hamiltonian, S_plus_1 , S_minus_N , S_Z_indexed matrices

H = sparse(0);

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
S_Z_indexed = cell(N,1);
S_number_indexed = cell(N,1);

for i=1:N
    S_Z_i = kron(speye(2^(i-1)),S_Z);
    S_Z_i = kron(S_Z_i,speye(2^(N-i)));
    S_number_i = kron(speye(2^(i-1)),S_plus*S_minus);
    S_number_i = kron(S_number_i,speye(2^(N-i)));
    S_Z_indexed{i} = S_Z_i;
    S_number_indexed{i} = S_number_i;
end

%% Initial state
State = [1;0];
for k = 2:N
    State = kron(State,[0;1]);
end
%%
Profiles_test = zeros(length(time),N);
for i = 1:length(time)
    
    for x = 1:N
        Profiles_test(i,x) = State'*((S_Z_indexed{x})*State);
    end
    
    %RK-4
    k1 = -(1i)*(H*State);
    k2 = -(1i)*H*(State + (dt/2)*k1);
    k3 = -(1i)*H*(State + (dt/2)*k2);
    k4 = -(1i)*H*(State + dt*k3);
    State = State + (dt/6)*(k1 + 2*k2 + 2*k3 + k4);
    State = State/norm(State);
    
end
%%
save('N = 9,Comparison','Profiles_test','N')
%% Animation
Profiles = Profiles_test;
x = 1:N;
y = Profiles(1,:);

for i = 2:length(Profiles)
y = Profiles(i,:);
plot(x,y)
axis([0.5 N+0.5 -1.2 1.2])
pause(1/30)
end

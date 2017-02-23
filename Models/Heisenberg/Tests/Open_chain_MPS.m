clc
close all
clear all
%%
N = 5;
J = 1;
U = 1;
G = 1;
dt = 0.01;
D_max = 15;
T = 1;

%rng shuffle
stream = RandStream('mt19937ar','Seed',5489); % MATLAB's start-up settings
RandStream.setGlobalStream(stream);
%% Stuff for conventional
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

State_O2 = random_mps(N,D_max,2);

[U_odd,U_even] = HeisenbergOpen_U_O2(N,J,U,G,dt);
%%
time = linspace(0,T,floor(T/dt)+1);
Profiles_MPS_O2 = zeros(length(time),N);
for i = 1 : length(time);
  
    
    %%%%%%%%%%%%%%%%%%%%%%%MPS%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%Order 2 Trotter%%%%%%%%%%%%%%%%
    
    for j = 1:N
        Sz_State = State_O2;
        Sz_State{j} = contract(State_O2{j},3,S_Z,2);
        Profiles_MPS_O2(i,j) = real(braket(Sz_State,State_O2));
    end
    
    S1 = State_O2;
    S1{1} = contract(State_O2{1},3,S_plus,2);
    P1 = real(G*dt*braket(S1,S1));
    
    S2 = State_O2;
    S2{N} = contract(State_O2{N},3,S_minus,2);
    P2 = real(G*dt*braket(S2,S2));
    
    %Stochastic evolution
    e = rand();
    
    if e <= P1
        State_O2 = sweep(S1,1);
    elseif e <= (P2 + P1) && e >= P1
        State_O2 = sweep(S2,1);
    else
        State_O2 = apply(U_odd,State_O2);
        State_O2 = sweep(State_O2,1);
        State_O2 = sweep(State_O2,-1,D_max);
        State_O2 =  apply(U_even,State_O2);
        State_O2 = sweep(State_O2,-1);
        State_O2 = sweep(State_O2,1,D_max);
        State_O2 = apply(U_odd,State_O2);
        State_O2 = sweep(State_O2,1);
        State_O2 = sweep(State_O2,-1,D_max);
    end
    
end
%%
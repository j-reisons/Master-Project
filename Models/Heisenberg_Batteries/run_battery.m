clc
close all
clear all
%%
N = 10;
J = 2;
U_c = 1;
U_b = 0;
d = 2;

D_start = 30;
ground_error = 1E-10;

T = 3;
dt = 0.05;
steps = round(T/dt);

D_max = 100;
comp_error = 1E-10;

tag = 'testing';

filename = ['Batteries','_N',strrep(num2str(N),'.',',') ,'_Ub',strrep(num2str(U_b),'.',',')...
    ,'_Uc',strrep(num2str(U_c),'.',','),'_','T',strrep(num2str(T),'.',','),'_','Dmax',num2str(D_max)...
    ,'_',tag,'.mat'];

%% Initial state preparation
State = cell(1,3*N);

Up = zeros(1,1,2);
Up(1,1,1) = 1;
Down = zeros(1,1,2);
Down(1,1,2) = 1;

H_Heis = Heisenberg_H(N,J,U_c);
MPS = random_mps(N,D_start,d,-1);
Ground = Iter_ground(H_Heis,MPS,ground_error);

for i=1:N
    State{i} = Up;
    State{N+i} = Ground{i};
    State{2*N + i} = Down;
end

%%
[U_even_dt,U_odd_dt] = Heisenberg_Batteries_U(N,J,U_b,U_c,dt);
[U_even_half,U_odd_half] = Heisenberg_Batteries_U(N,J,U_b,U_c,dt/2);
U = compressMPO(U_odd_half,U_even_dt,U_odd_half);

Magnetizations = zeros(3*N,steps+1);
Currents = zeros(3*N,steps);

S_Z =[
    [1 , 0]
    [0 ,-1]
    ];

%%
fprintf('Three U \n')
tic
for j = 1:3*N
    Sz_State = State;
    Sz_State{j} = contract(State{j},3,S_Z,2);
    Magnetizations(j,1) = real(braket(Sz_State,State));
end

for i = 1:steps
    
    for j = 1:3*N
        Sz_State = State;
        Sz_State{j} = contract(State{j},3,S_Z,2);
        Magnetizations(j,i+1) = real(braket(Sz_State,State));
    end
    
    State = apply(U_odd_half,State);
    State = Iter_comp(State,D_max,comp_error);
    
    State = apply(U_even_dt,State);
    State = Iter_comp(State,D_max,comp_error);
    
    State = apply(U_odd_half,State);
    State = Iter_comp(State,D_max,comp_error);
end
toc 
%%
save(filename,'Magnetizations')
%%

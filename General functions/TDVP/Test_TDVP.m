clc
close all
clear all
%%
% I run Heisenberg Batteries model for known parameters and compare
% currents and magnetization results between TEBD and TDVP

N = 10;
Uc = 0.5;
Ub = 0;
tag = 'TDVP';

J = 1;
d = 2;

D_start = 100;
D_max = 200;

T = N/(2*J);
dt = 0.05/J;
steps = round(T/dt);

ground_error = 1E-8;

filename = ['Batteries','_N',strrep(num2str(N),'.',',') ,'_Ub',strrep(num2str(Ub),'.',',')...
    ,'_Uc',strrep(num2str(Uc),'.',','),'_','T',strrep(num2str(T),'.',','),'_','Dmax',num2str(D_max)...
    ,'_',tag,'.mat'];

%% Initial state preparation

State = cell(1,3*N);

Up = zeros(1,1,2);
Up(1,1,1) = 1;
Down = zeros(1,1,2);
Down(1,1,2) = 1;

H_Heis = Heisenberg_H(N,J,Uc);
MPS = random_mps(N,D_start,d,-1);
Ground = Iter_ground(H_Heis,MPS,ground_error);

for i=1:N
    State{i} = Up;
    State{N+i} = Ground{i};
    State{2*N + i} = Down;
end

%%
Hamiltonian = Batteries_H(N,J,Uc,Ub);

Magnetizations = zeros(3*N,steps+1);
Currents = zeros(3*N - 1,steps+1);

S_Z = [1,0;0,-1];
S_X =[0,1;1,0];
S_Y = [0,-1i;1i,0];

S_Z_mpo = cell(1,1);
S_Z_1 = reshape(S_Z,[1,1,d,d]);
S_Z_mpo{1} = S_Z_1;

Q_mpo = cell(1,2);  
Q_1 = zeros(1,2,d,d);
Q_1(1,1,:,:) = S_X*((2*J)^0.5);
Q_1(1,2,:,:) = S_Y*((2*J)^0.5);
Q_mpo{1} = Q_1;
Q_2 = zeros(2,1,d,d);
Q_2(1,1,:,:) = S_Y*((2*J)^0.5);
Q_2(2,1,:,:) = - S_X*((2*J)^0.5);
Q_mpo{2} = Q_2;

%%
State = Inflate_mps(State,D_max);
State = sweep(State,-1);
canon = -1;
%%

for i = 1:steps
    evaluations = Canon_evaluator(State,canon,S_Z_mpo,Q_mpo);
    Magnetizations(:,i) = real(evaluations{1});
    Currents(:,i) = real(evaluations{2});
    
    State = Iter_evolve(State,Hamiltonian,dt);
    State = sweep(State,-1);
end

evaluations = Canon_evaluator(State,canon,S_Z_mpo,Q_mpo);
Magnetizations(:,steps+1) = real(evaluations{1});
Currents(:,steps+1) = real(evaluations{2});

%%
save(filename,'Magnetizations','Currents','State')

%%
load('Batteries_N10_Ub0_Uc0,5_T5_Dmax100_TDVP.mat');
Magnetizations_TDVP = Magnetizations;
Currents_TDVP = Currents;
load('Batteries_N10_Ub0_Uc0,5_T5_Dmax100_threeU.mat')

%%
imagesc(Currents-Currents_TDVP)
figure
imagesc(Magnetizations-Magnetizations_TDVP)
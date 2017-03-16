clc
close all
clear all
%%
N = 10;
Uc = 0.5;
Ub = 0;

J = 1;
d = 2;

D_start = 100;
D_max = 100;

ground_error = 1E-8;

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
fprintf('Testing Inflation...')
FatState = Inflate_mps(State,200);
FatState = sweep(FatState,1);
FatState = sweep(FatState,-1);

assert(approx(braket(FatState,State),1,1E-14));

fprintf('done.\n')
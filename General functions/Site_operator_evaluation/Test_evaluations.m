clc
close all
clear all
%%
N = 30;
J = 1;
U = 1;
D_start = 40;
d = 2;
ground_error = 1E-7;

%%
State = cell(1,3*N);

Up = zeros(1,1,2);
Up(1,1,1) = 1;
Down = zeros(1,1,2);
Down(1,1,2) = 1;

H_Heis = Heisenberg_H(N,J,U);
MPS = random_mps(N,D_start,d,-1);
Ground = Iter_ground(H_Heis,MPS,ground_error);

for i=1:N
    State{i} = Up;
    State{N+i} = Ground{i};
    State{2*N + i} = Down;
end


%%
S_Z = [1,0;0,-1];
S_X =[0,1;1,0];
S_Y = [0,-1i;1i,0];

S_Z_mpo = cell(1,1);
S_Z_1 = reshape(S_Z,[1,1,d,d]);
S_Z_mpo{1} = S_Z_1;

ZZ_mpo = cell(1,2);
ZZ_mpo{1} = reshape(S_Z,[1,1,d,d]);
ZZ_mpo{2} = ZZ_mpo{1};
%%
Magnet_braket = zeros(1,3*N);
Magnet_contract = zeros(1,3*N);
Magnet_mixed_can = zeros(1,3*N);


%%

fprintf('Testing braket evaluation \n')
tic
for i = 1:3*N
    Sz_State = State;
    Sz_State{i} = contract(State{i},3,S_Z,2);
    Magnet_braket(i) = real(braket(Sz_State,State));
end
toc
   

State_left = sweep(State,1);
fprintf('Testing Iterative evaluation \n')
tic
R = R_contractions(State_left);
for i = 1:3*N
    Magnet_L_site = contract(R{i},1,State_left{i},2);
    Magnet_L_site = contract(Magnet_L_site,3,S_Z,2);
    Magnet_L_site = contract(Magnet_L_site,[1,3],conj(State_left{i}),[2,3]);
    Magnet_contract(i) = trace(Magnet_L_site);
end
toc

fprintf('Testing mixed-canonized evaluation (Right) \n')
State_right = sweep(State,-1);

tic
Magnet_canon = Canon_evaluator(State_right,-1,S_Z_mpo,ZZ_mpo);
Magnet_canon_right = Magnet_canon{1};
Z_Z_right = Magnet_canon{2};
toc

fprintf('Testing mixed-canonized evaluation (Left) \n')
State_left = sweep(State,1);

tic
Magnet_canon = Canon_evaluator(State_left,1,S_Z_mpo,ZZ_mpo);
Magnet_canon_left = Magnet_canon{1};
current_canon_left = Magnet_canon{2};
toc

%%

assert(approx(max(Magnet_braket - Magnet_contract),0,1E-10));
assert(approx(max(Magnet_braket - Magnet_canon_left),0,1E-10));
assert(approx(max(Magnet_braket - Magnet_canon_right),0,1E-10));




close all;
clear all;
clc

N = 10;     % Number of sites
d = 2;      % Local H-space dimension
D_S = 2;    % State bond dim
D_O = 3;    % Operator bond dim

%% MPS initialization
M = random_mps(N,10,2);
%% MPO initialization

S_z = [1,0;0,-1];
W = zeros(D_O,D_O,d,d);

W(1,1,:,:) = eye(2);
W(2,1,:,:) = S_z;
W(3,2,:,:) = S_z;
W(3,3,:,:) = eye(2);

O = cell(1,N);
O{1} = W(D_O,:,:,:);
for i = 2:N-1
    O{i} = W;
end
O{N} = W(:,1,:,:);

%% Testing canonization/compression functions
Mprime = apply(O,M);
newstate = expand_MPS(Mprime);

Mleft = Mprime;
Mright = Mprime;
for i = 1:N
Mleft = L_can(Mleft,i);
end
for i = N:-1:1
Mright = R_can(Mright,i);
end

expand_MPS(Mleft)'*(newstate/norm(newstate)) % + + + + + + + + + + + + + +
%Mright seems to be sometimes +, sometimes - the state,while Mleft is
%always + the state
expand_MPS(Mright)'*(newstate/norm(newstate)) % - - + - - + - - + + - - + +

%% Testing compression
M_r_to_l = Mright;
M_l_to_r = Mleft;

for i = 1:N
M_r_to_l = L_can(M_r_to_l,i,3);
end

for i = N:-1:1
M_l_to_r = R_can(M_l_to_r,i,3);
end

expand_MPS(M_r_to_l)'*(newstate/norm(newstate))

expand_MPS(M_l_to_r)'*(newstate/norm(newstate))
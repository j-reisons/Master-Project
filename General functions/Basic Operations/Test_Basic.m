close all;
clear all;

N = 10;     % Number of sites
d = 2;      % Local H-space dimension
D_S = 2;    % State bond dim
D_O = 3;    % Operator bond dim

%% MPS initialization
M = cell(1,N);

A_odd = zeros(D_S,D_S,d);   A_odd(:,:,1) = eye(D_S);
A_even = zeros(D_S,D_S,d);  A_even(:,:,2) = eye(D_S);

M{1} = zeros(1,2,2);
M{1}(:,:,1) = [1,0];    M{1}(:,:,2) = [0,0];

for i = 2:2:N-1
    M{i} = A_even;
    M{i + 1} = A_odd;
end

M{N} = zeros(2,1,2);
if mod(N,2) == 1
    M{N}(:,:,1) = [1;0];    M{N}(:,:,2) = [0;0];
else
     M{N}(:,:,1) = [0;0];    M{N}(:,:,2) = [1;0];
end
    
%% MPO initialization
S_z = 1/sqrt(2)*[1,0;0,-1];
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

%% Traditional matrix operations
neel = [1;0];
for i = 2:N
    if mod(i,2) == 0;
        neel = kron(neel,[0;1]);
    else
        neel = kron(neel,[1;0]);
    end
end

ham = zeros(2^N);
for i = 1:(N-1) % number of terms in H
    S = sparse(1);
    for j = 1:N 
        if j == i
            S = kron(S,S_z);
        elseif j == i + 1
            S = kron(S,S_z);
        else
            S = kron(S,speye(2));
        end
    end
    ham = ham + S;
end

%% Calculate new state
Mprime = apply(O,M);
newstate = expand_MPS(Mprime);
%%
fprintf('Testing MPO*MPS...');
assert(isequal(newstate,ham*neel))
fprintf('\t\tdone\n');
%%
fprintf('Testing scalar product...');
assert(isequal(braket(M,Mprime),neel'*ham*neel))
fprintf('\tdone\n');
%%
bigO = expand_MPO(O);
tolerance = 1e-15;

fprintf('Testing MPO expansion...');
assert(approx(max(max(bigO-ham)),0,tolerance));
fprintf('\tdone\n');

%% Testing canonized state

Mleft = Mprime;
Mleft = sweep(Mprime,1);


overlap = expand_MPS(Mleft)'*(newstate/norm(newstate));

fprintf('Testing left canonization...');
assert(approx(overlap,1,tolerance))
fprintf('\tdone\n');

function [H] = Batteries_H (N,J,Uc,Ub)
% Batteries_H Returns MPO representation of XXZ B
% Procedure outline in Schollwock 6.1

H = cell(1,3*N);

S_X = [0,1;1,0];
S_Y = [0,-1i;1i,0];
S_Z = [1,0;0,-1];
I = eye(2);


left = zeros(1,5,2,2);
bulk_chain = zeros(5,5,2,2);
bulk_bath = zeros(5,5,2,2);
right = zeros(5,1,2,2);

left(1,2,:,:) = J*S_X;
left(1,3,:,:) = J*S_Y;
left(1,4,:,:) = Ub*S_Z;
left(1,5,:,:) = I;

bulk_chain(1,1,:,:) = I;
bulk_chain(2,1,:,:) = S_X;
bulk_chain(3,1,:,:) = S_Y;
bulk_chain(4,1,:,:) = S_Z;

bulk_chain(5,2,:,:) = J*S_X;
bulk_chain(5,3,:,:) = J*S_Y;
bulk_chain(5,4,:,:) = Uc*S_Z;
bulk_chain(5,5,:,:) = I;

bulk_bath = bulk_chain;
bulk_bath(5,4,:,:) = Ub*S_Z;

right(1,1,:,:) = I;
right(2,1,:,:) = S_X;
right(3,1,:,:) = S_Y;
right(4,1,:,:) = S_Z;

H{1} = left;
for i = 2:N
    H{i} = bulk_bath;
end

for i = N+1:(2*N)-1
    H{i} = bulk_chain;
end

for i = 2*N:(3*N)-1
    H{i} = bulk_bath;
end

H{3*N} = right;
end


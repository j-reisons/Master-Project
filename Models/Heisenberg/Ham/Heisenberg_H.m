function [H] = Heisenberg_H (N,J,U)
% Heisenberg_H Returns MPO representation of XXZ hamiltonian
% Procedure outline in Schollwock 6.1

H = cell(1,N);

S_X = [0,1;1,0];
S_Y = [0,-1i;1i,0];
S_Z = [1,0;0,-1];
I = eye(2);


left = zeros(1,5,2,2);
bulk = zeros(5,5,2,2);
right = zeros(5,1,2,2);

left(1,2,:,:) = J*S_X;
left(1,3,:,:) = J*S_Y;
left(1,4,:,:) = U*S_Z;
left(1,5,:,:) = I;

bulk(1,1,:,:) = I;
bulk(2,1,:,:) = S_X;
bulk(3,1,:,:) = S_Y;
bulk(4,1,:,:) = S_Z;

bulk(5,2,:,:) = J*S_X;
bulk(5,3,:,:) = J*S_Y;
bulk(5,4,:,:) = U*S_Z;
bulk(5,5,:,:) = I;

right(1,1,:,:) = I;
right(2,1,:,:) = S_X;
right(3,1,:,:) = S_Y;
right(4,1,:,:) = S_Z;

H{1} = left;
for i = 2:N-1
    H{i} = bulk;
end
H{N} = right;

end


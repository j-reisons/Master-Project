function operator = expand_MPO(mpo)
%expands MPO into "Big matrix"
%   Detailed explanation goes here

N = length(mpo);
d = size(mpo{1},4);

middle = floor(N/2);
parity = mod(N,2);

L = contract(mpo{2},1,mpo{1},2);
for i = 3: middle
    L = contract(mpo{i},1,L,1);
end
L = squeeze(L);

R = contract(mpo{N},1,mpo{N-1},2);
for i = N - 2: -1 : middle + 1
    R = contract(R,2*(N-i),mpo{i},2);
end
R = squeeze(R);

sR = size(R);

operator = contract(R,length(sR)-2,L,1);

operator = permute(operator,[1:2:2*N,2:2:2*N]);
operator = reshape(operator,[d^N,d^N]);
end


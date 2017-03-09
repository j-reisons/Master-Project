function operator = expand_MPO(mpo)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

N = length(mpo);
d = size(mpo{1},4);

middle = floor(N/2);

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
sL = size(L);
if length(sR) > 2*(middle+1)
operator = contract(R,length(sR)-2,L,1);
else
    R = reshape(R,[1,sR]);
    L = reshape(L,[1,sL]);
    operator = contract(R,1,L,1);
end
operator = permute(operator,[1:2:2*N,2:2:2*N]);
operator = reshape(operator,[d^N,d^N]);
end


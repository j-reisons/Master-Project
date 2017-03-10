function [state] = expand_MPS(mps)
% expands MPS into "Big vector". Does not work for N < 4
%   Detailed explanation goes here

N = length(mps);
d = size(mps{1},3);
middle = floor(N/2);

L = contract(mps{2},1,mps{1},2);
if middle > 3
    for i = 3: middle
        L = contract(mps{i},1,L,1);
    end
end

L = squeeze(L);

R = contract(mps{N},1,mps{N-1},2);
R = squeeze(R);


for i = N - 2: -1 : middle + 1
    R = contract(R,N-i,mps{i},2);
end


R = squeeze(R);

sR = size(R);
state = contract(R,length(sR)-1,L,1);
state = reshape(state,[d^N,1]);

end


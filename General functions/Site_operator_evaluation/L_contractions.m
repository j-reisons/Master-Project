function [L] = L_contractions(mps)
%L_CONTRACTIONS Summary of this function goes here
%   Detailed explanation goes here

N = length(mps);
L = cell(1,N);
L{1} = 1;

for i = 2:1:N
    L{i} = contract(L{i-1},1,mps{i-1},1);
    L{i} = contract(L{i},[1,3],conj(mps{i-1}),[1,3]);
end

end


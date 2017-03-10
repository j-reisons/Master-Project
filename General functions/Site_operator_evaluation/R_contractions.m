function [R] = R_contractions(mps)
%R_CONTRACTIONS Summary of this function goes here
%   Detailed explanation goes here

N = length(mps);
R = cell(1,N);
R{N} = 1;

for i = N-1:-1:1
    R{i} = contract(R{i+1},1,mps{i+1},2);
    R{i} = contract(R{i},[3,1],conj(mps{i+1}),[3,2]);
end
end


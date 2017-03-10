function [sprod] = worse_braket(mps_1,mps_2)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

N = length(mps_1);
d = size(mps_1{1},3);

contraction = 1;

for i = 1:N
    contraction = contract(contraction,1,mps_2{i},1);
    contraction = contract(contraction,[1,3],conj(mps_1{i}),[1,3]);
end
contraction = trace(contraction);

sprod = contraction;


end


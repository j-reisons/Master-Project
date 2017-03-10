function sprod = braket(mps_1,mps_2)
% Returns the scalar product < MPS_1 | MPS_2 >.
% The inputs are cell arrays corresponding to MPS decompositions
% I don't quite remember writing this, but IIRC it follows the procedure of
% Scholl 4.2.1 page 34

N = length(mps_1);
d = size(mps_1{1},3); % Physical dimension
sprod = 0;

for sigma = 1:d % First site
    sprod = sprod + mps_1{1}(:,:,sigma)'*mps_2{1}(:,:,sigma);
end

for i = 2:N
    sum = 0;
    for sigma= 1:d
        sum = sum + mps_1{i}(:,:,sigma)'*sprod*mps_2{i}(:,:,sigma);
    end
    sprod = sum;
end

end
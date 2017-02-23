function state = expand_MPS(mps)
% expands MPS into "big vector"
% Giacommo wrote this version of the function, I can't comment on it

N = length(mps);
d = size(mps{1},3);

state = zeros(d^N,1);
c = cell(1,N);  % The following code expands the MPS, may take a lot of mem
[c{:}] = ndgrid(1:d);
combs = fliplr(cell2mat(cellfun(@(v)v(:),c,'UniformOutput',false)));
for pos = 1:d^N
    prod = 1;
    for i = 1:N
        prod = prod*mps{i}(:,:,combs(pos,i));
    end
    state(pos) = prod;
end

end


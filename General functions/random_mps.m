function [mps] = random_mps(N,D,d,varargin)
%RANDOM_MPS random mps of spin chain length N, bond length D, local H space
%dimension d. Returns Right or left normalized if specified, by default
%left

mps = cell(1,N);
mps{1} = rand(1,D,d);
for i = 2:N-1
   mps{i} = rand(D,D,d);
end
mps{N} = rand(D,1,d);

direction = 1;
if ~isempty(varargin)
    direction = varargin{1};
end

switch direction
    case 1
        mps = sweep(mps,direction);
        return
    case -1
        mps = sweep(mps,direction);
        return
end
end


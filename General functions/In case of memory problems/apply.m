function mps = apply(mpo,mps,varargin)
% Applies an MPO to an MPS. Can apply mps on only one site by setting
% varargin = site. If unspecified, applies on the whole chain
%
% The MPO is assumed to be the same for each site
% i.e. a D_O*D_O square matrix of operators 
% where D_O is the operator bond dimension.
% order of indices is (bond,bond,physical,physical)
% Boundary conditions are assumed to be open.
%
% The MPS is assumed to be a cell array
%
% See Schollwock 5.1 p.56 

if ~isempty(varargin)
    site = varargin{1};
    result = contract(mps{site},3,mpo{site},4);
    R_size = size(result);
    result = reshape(permute(result,[3 1 4 2 5]),[R_size(1)*R_size(3),R_size(2)*R_size(4),R_size(5)]);
    mps{site} = result;
    return
end

N = length(mps);
result = cell(1,N);

for i = 1:N
    result{i} = contract(mps{i},3,mpo{i},4);
end

for i = 1:N 
        R_size = size(result{i});
        result{i} = reshape(permute(result{i},[3 1 4 2 5]),[R_size(1)*R_size(3),R_size(2)*R_size(4),R_size(5)]);
end
mps = result;
end
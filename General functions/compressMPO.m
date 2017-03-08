function [U,U_norm] = compressMPO(varargin)
% Please don't put non-MPOS into varargin
% Beware of ordering of MPOS as arguments, as they don't commute

for n_mpo = 1:length(varargin)
	if ~iscell(varargin{n_mpo})
		break;
	end
end

N = length(varargin{1});
D = size(varargin{1}{1},3);

U = cell(1,N);
for site = 1:N
		W = varargin{1}{site};
	if n_mpo > 1
		for k = 2:n_mpo
			W = multiply(W,varargin{k}{site});
		end
	end
    [b_1,b_2,~,~] = size(W);
    U{site} = reshape(W,[b_1,b_2,D*D]);
end

[U,U_norm] = sweep(U,-1);
[U,U_sign] = sweep(U,1,eps,1000);
U_norm = real(U_norm*U_sign);

for site = 1:N
	U{site} = U_norm^(1/N)*U{site};
end

for site = 1:N
    [b_1,b_2,~,~] = size(U{site});
    U{site} = reshape(U{site},[b_1 b_2 D D]);
end
end
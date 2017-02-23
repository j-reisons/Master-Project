function [Mleft,error,D_max] = L_can(mps,site,varargin)
% Left canonizes (if asked compresses) $site of mps and throws the SV' to the right
% varargin is empty, or tolerance and D_limit
% If empty, canonization without compression
% If tolerance and D_lmit is provided compresses to tolerance and D_limit
% Last site normalization is thrown out
% Procedure is outline in Schollwock 4.4.1 , page 43

a = 0; % 0 is no compression,1 is compression

if ~isempty(varargin)
    tolerance = varargin{1};
    D_limit = varargin{2};
    a = 1; % a = 1 is compression
else
    a = 0 ; % a = 0 is no compression
end

N = length(mps);
Mleft = mps; % returned MPS

% Middle sites    
work = Mleft{site}; % Intermediate "work" variable
s_w = size(work); % Will be needed later for a reshape
% Reshaping to perform SVD
work = permute(work,[1 3 2]);
work = reshape(work,[s_w(3)*s_w(1),s_w(2)]);

try
    [U,S,V] = svd(work,'econ');
catch % SVD not-converging workaround
    work = work + rand(size(work))*1E-10;
    [U,S,V] = svd(work,'econ');
end
s_s = size(S,1); % Schmidt rank

switch a
    case 0
    S = S /sqrt(trace(S*S')); % Imposing sum of squares = 1
    D_max = size(S,1);        % Return new bond dimension
    error = 0;                % error = 0 ( no compression)
    case 1
    S_2 = S*S';
    S = S /sqrt(trace(S_2)); % Imposing sum of squares = 1
    S_2 = diag(S*S');
    cut = 0;
    sum = 0;
    % We find the dimension to which we need to compress to respect tolerance
    while 1 - sum > tolerance && cut < D_limit && cut < length(S_2)
        cut = cut +1;
        sum = sum + S_2(cut);
    end
    error = 1 - sum;
    D_max = cut;
    % Compression
    U = U(:,1:D_max);
    S = S(1:D_max,1:D_max);
    S = S /sqrt(trace(S*S')); % Imposing sum of squares = 1
    V = V(:,1:D_max);
end
s_u = size(U);
A = reshape(U,[s_w(1) s_w(3) s_u(2)]);
A = permute(A,[1 3 2]);
Mleft{site} = A;
SV = S*V';
if site ~= N
    Mleft{site+1} = contract(SV,2,Mleft{site+1},1);
else
    Mleft{site} = Mleft{site}*sign(SV); % Fixing +- 1 phase (I noticed I would
    %sometimes get the state with a wrong -1 phase. I don't know why, but
    %this fixes it)
end

end
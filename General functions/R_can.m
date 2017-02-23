function [Mright,error,D_max] = R_can(mps,site,varargin)
% Right canonizes (if asked compresses) $site of mps and throws the US to the left.
% varargin is empty, or tolerance and D_limit
% If empty, canonization without compression
% If tolerance and D_lmit is provided compresses to tolerance and D_limit
% First site normalization is thrown out
% Procedure is outline in Schollwock 4.4.2 , page 44

if ~isempty(varargin)
    tolerance = varargin{1};
    D_limit = varargin{2};
    a = 1; % a = 1 is compression
else
    a = 0; % a = 0 is no compression
end

N = length(mps);
Mright = mps; % returned MPS

work = Mright{site}; % Intermediate "work" variable
s_w = size(work);
work = reshape(work,[s_w(1),s_w(2)*s_w(3)]);
try
[U,S,V] = svd(work,'econ');
catch % SVD not-converging workaround
    work = work + rand(size(work))*1E-10;
    [U,S,V] = svd(work,'econ');
end
s_s = size(S,1);
switch a
    case 0
        S = S /sqrt(trace(S*S'));
        D_max = size(S,1);
        error = 0;
    case 1
    S_2 = S*S';
    S = S /sqrt(trace(S_2));
    S_2 = diag(S*S');
    cut = 0;
    sum = 0;
    while 1 - sum > tolerance && cut < D_limit && cut < length(S_2)
        cut = cut +1;
        sum = sum + S_2(cut);
    end
    error = 1 - sum;
    D_max = cut;
    U = U(:,1:D_max);
    S = S(1:D_max,1:D_max);
    S = S /sqrt(trace(S*S'));
    V = V(:,1:D_max);
end
s_v = size(V');
B = reshape(V',[s_v(1),s_w(2),s_w(3)]);
Mright{site} = B;
US = U*S;
if site ~=1
    Mright{site-1} = contract(Mright{site-1},2,US,1);
    Mright{site-1} = permute(Mright{site-1},[1,3,2]);
else
    Mright{site} = Mright{site}*sign(US); % Fixing +- 1 phase
end

end


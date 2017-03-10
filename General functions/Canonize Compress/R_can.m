function [Mright,mps_norm] = R_can(mps,site,varargin)
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

mps_norm = 0;
N = length(mps);
Mright = mps; % returned MPS

work = Mright{site}; % Intermediate "work" variable
s_w = size(work);
work = reshape(work,[s_w(1),s_w(2)*s_w(3)]);

switch a
    case 0
        
        [Q,R] = qr(work',0);
        US = R';
        V_dag = Q';
        
    case 1
        
        try
            [U,S,V] = svd(work,'econ');
        catch % SVD not-converging workaround
            work = work + rand(size(work))*1E-10;
            [U,S,V] = svd(work,'econ');
        end
        
        S_2 = diag(S*S');
        S_2 = S_2/sum(S_2);
        
        cut = 0;
        running_sum = 0;
        % We find the dimension to which we need to compress to respect tolerance
        while 1 - running_sum > tolerance && cut < D_limit && cut < length(S_2)
            cut = cut + 1;
            running_sum = running_sum + S_2(cut);
        end
        
        % Compression
        U = U(:,1:cut);
        S = S(1:cut,1:cut);
        V = V(:,1:cut);
        
        US = U*S;
        V_dag = V';
end

s_v = size(V_dag);
B = reshape(V_dag,[s_v(1),s_w(2),s_w(3)]);
Mright{site} = B;

if site ~=1
    Mright{site-1} = contract(Mright{site-1},2,US,1);
    Mright{site-1} = permute(Mright{site-1},[1,3,2]);
else
    Mright{site} = Mright{site}*sign(US);     % Fixing +- 1 phase
    mps_norm = sign(US)*US;

end

end


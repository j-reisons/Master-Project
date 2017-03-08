function [mps_out,norm] = sweep(mps_in,direction,varargin)
% Function canonizes (if asked compresses) entire state (it sweeps through
% it)
% varargin is empty, or tolerance and D_limit
% If empty, canonization without compression
% If tolerance and D_lmit is provided compresses to tolerance or D_limit
% tolerance/D_limit are simply passed along to R_can and L_can, 
% nothing happens in this function

mps_out = mps_in;
N = length(mps_in);

switch direction
    case 1
        for i = 1:N
            [mps_out,norm] = L_can(mps_out,i,varargin{:});
        end
        
    case -1
        for i = N:-1:1
            [mps_out,norm] = R_can(mps_out,i,varargin{:});
        end
        
end
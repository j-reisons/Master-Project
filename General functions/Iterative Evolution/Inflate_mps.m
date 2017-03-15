function [fat_mps] = Inflate_mps(mps,D_max)
%INFLATE_MPS Artificially inflates mps to D_max size, for a dirty starting
% point to 1-site iterative evolution
mps = sweep(mps,1);
mps = sweep(mps,-1);

N = length(mps);
s = size(mps{1});
d = s(3);
fat_mps = cell(1,N);
N_growth = floor(log(D_max)/log(d));
smalltrash = 1E-12;

if N_growth >= floor(N/2)
    a = 1;
else
    a = 2;
end

switch a
    case 1
        
        for i = 1:floor(N/2)
            M = rand(d^(i-1),d^i,d)*smalltrash;
            s = size(mps{i});
            M(1:s(1),1:s(2),1:s(3)) = mps{i};
            fat_mps{i} = M;
            
            M = rand(d^i,d^(i-1),d)*smalltrash;
            s = size(mps{N-i+1});
            M(1:s(1),1:s(2),1:s(3)) = mps{N-i+1};
            fat_mps{N-i+1} = M;
        end
        if mod(N,2)
            M = rand(d^(floor(N/2)),d^(floor(N/2)),d)*smalltrash;
            s = size(mps{ceil(N/2)});
            M(1:s(1),1:s(2),1:s(3)) = mps{ceil(N/2)};
            fat_mps{ceil(N/2)} = M;
        end
        
    case 2
        
        for i = 1:N_growth
            M = rand(d^(i-1),d^i,d)*smalltrash;
            s = size(mps{i});
            M(1:s(1),1:s(2),1:s(3)) = mps{i};
            fat_mps{i} = M;
            
            M = rand(d^i,d^(i-1),d)*smalltrash;
            s = size(mps{N-i+1});
            M(1:s(1),1:s(2),1:s(3)) = mps{N-i+1};
            fat_mps{N-i+1} = M;
        end
        
        M = rand(d^N_growth,D_max,d)*smalltrash;
        s = size(mps{N_growth+1});
        M(1:s(1),1:s(2),1:s(3)) = mps{N_growth+1};
        fat_mps{N_growth + 1} = M;
        
        M = rand(D_max,d^N_growth,d)*smalltrash;
        s = size(mps{N - N_growth});
        M(1:s(1),1:s(2),1:s(3)) = mps{N - N_growth};
        fat_mps{N - N_growth} = M;
        
        for i = N_growth+2 : N - N_growth -1
            M = rand(D_max,D_max,d)*smalltrash;
            s = size(mps{i});
            M(1:s(1),1:s(2),1:s(3)) = mps{i};
            fat_mps{i} = M;
        end
        
end

end


function [ft_X , Omega] = dumbFT(X,dt,dw)
%DUMBFT Summary of this function goes here
%   Detailed explanation goes here
    N = length(X);
    Omega = 0:dw:dw*(N-1);

    M = zeros(N,N);
    
    for m = 1:N
        for n = 1:N
            M(m,n) = exp(-1i*dt*dw*(m-1)*(n-1));
        end
    end
    
    ft_X = M*(X.');
end


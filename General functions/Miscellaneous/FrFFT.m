function [ft_X , Omega] = FrFFT(X,dt,dw)
% Returns

    gammaFrFFT = dt*dw/(2*pi);
    N = length(X);
    Omega = [0:dw:dw*(N-1)];
    
    y = zeros(2*N,1);
    y(1:N) = X(:).*exp(-1i*pi*gammaFrFFT*(((1:N)-1).^2)).';
    
    z = zeros(2*N,1);
    z(1:N) = exp(1i*pi*gammaFrFFT*(((1:N)-1).^2));
    z(N:2*N) = exp(1i*pi*gammaFrFFT*((2*N-(N:2*N)+1).^2));
    
    oversized = exp(-1i*pi*gammaFrFFT*( ((1:2*N)-1).^2 ) )...
        .*ifft(fft(z).*fft(y)).';
    ft_X = oversized(1:N);
    
end
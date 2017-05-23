clc
close all
clear all
%%
N = 20;
N = 4*N;
D = 100;
d = 2;

%Tits

U_multiplier = 4*16;
complex_number = 16;
%%
n = floor(log2(D));
if 2*n >= N;
    n = floor(N/2);
    m = 0;
else
    m = N - 2*n; 
end
%%
Big = m*d*D*D;
Tails = 0;
for i = 1:n
    Tails = Tails + 2*d^(1+2*(i-1));
end
MPS = complex_number*(Tails + Big);

%% Iterative compression
L_R = 2*(10/d)*MPS;
out = MPS;
comp = L_R + out;

%% Total
Total = 64*MPS + comp;
Cheap = MPS + comp;

%%
fprintf(['This MPS will take ',num2str(MPS/1E9),' GB before evolution, and ',num2str(64*MPS/1E9),...
    ' after evolution. \n']);
fprintf(['Expensive evolution should take ' , num2str(Total/1E9), ' GB \n']);
fprintf(['Cheap evolution should take ' , num2str(Cheap/1E9), ' GB \n']);
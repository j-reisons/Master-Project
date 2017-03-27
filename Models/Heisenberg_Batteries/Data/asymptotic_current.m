clc
close all
clear all
%%
J= 1;
U = 0.5:0.1:1.5;

N = 20:10:50;
Ncs = N;
Nbs = 1.5*Ncs;
dt = 0.05;

D_max = 100;
freesweeps = 2;
cut = 0.8;

tag = '';
%% Filenames et al.

col=hsv(2*length(N));
asymptotics = zeros(1,length(U));
stdevs = zeros(1,length(U));

for i = 1:length(U)
    asymptotic = 0;
    squares = 0;
    for j = 1:length(N);  
        filename = ['Batteries','_Nc',strrep(num2str(Ncs(j)),'.',','),'_Nb',...
            strrep(num2str(Nbs(j)),'.',','),'_Ub',strrep(num2str(U(i)),'.',',')...
            ,'_Uc',strrep(num2str(U(i)),'.',','),'_','dt',strrep(num2str(dt),'.',',')...
            ,'_Dmax',num2str(D_max),'_f',num2str(freesweeps),'_',tag,'.mat'];
        
        load(filename)
        s = size(Currents);
        asymptotic = asymptotic + Currents(Nbs(j),round(cut*s(2)));
        squares = squares + Currents(Nbs(j),round(cut*s(2)))^2;
    end
    asymptotics(i) = asymptotic/length(N);
    stdevs(i) = sqrt(squares/length(N) - asymptotics(i)^2);

end

%% Charge plot
errorbar(U,asymptotics,stdevs);


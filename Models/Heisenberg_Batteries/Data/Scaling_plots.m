clc
close all
clear all
%%
J= 1;
tag = 'oneU';
U_b = 0;
U_c = 1.1;
D_max = 400;

%% Filename

Charges = cell(1);
Interfaces = cell(1);
Times = cell(1);

i = 1;
for N = 20:10:40
    
T = N/(2*J);
    
filename = ['Batteries','_N',strrep(num2str(N),'.',',') ,'_Ub',strrep(num2str(U_b),'.',',')...
    ,'_Uc',strrep(num2str(U_c),'.',','),'_','T',strrep(num2str(T),'.',','),'_','Dmax',num2str(D_max)...
    ,'_',tag,'.mat'];

load(filename);

% Magnetization charge
Charge = sum(Magnetizations(1:N,:))/N;
Charges{i} = Charge;

% Interface current
Interface = Currents(N,:);
Interfaces{i} = Interface;

% Rescaled times
Time = linspace(0,0.5,(10*N)+1);
Times{i} = Time;

i = i+1;
end

%% Charge plot
figure
for i = 1:3
    plot(Times{i},Charges{i})
    hold on
end

%% Current plot
figure
for i = 1:3
    plot(Times{i},Interfaces{i})
    hold on
end


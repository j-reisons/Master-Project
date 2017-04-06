clc
close all
clear all
%%
J= 1.2;
U = 1.2;

N = 50;
Ncs = N;
Nbs = 150%1.5*Ncs;

Cut = 4*Ncs;
Steps = Nbs*10-Cut;

dt = 0.05;
dw = 12/Steps;


TotalTime = 0:dt:10*Nbs*dt;
Time = 0:dt:Steps*dt;
Omega = 0:dw:Steps*dw;
Frequencies = Omega/(2*pi);

D_max = 100;
freesweeps = 2;
cut = 0.8;

tag = '';
%% Filenames et al.

col=hsv(length(U));
Fouriers = cell(1,1);
Currents_middle = cell(1,1);
Treated_middle = cell(1,1);
dumb_Fourier_middle = cell(1,1);
FRFFT_middle = cell(1,1);

for i = 1:length(U);
    
filename = ['Batteries','_Nc',strrep(num2str(Ncs),'.',','),'_Nb',...
    strrep(num2str(Nbs),'.',','),'_Ub',strrep(num2str(U(i)),'.',',')...
    ,'_Uc',strrep(num2str(U(i)),'.',','),'_','dt',strrep(num2str(dt),'.',',')...
    ,'_Dmax',num2str(D_max),'_f',num2str(freesweeps),'_',tag,'.mat'];

load(filename);

Current_middle = Currents(round(Nbs + (Ncs/2)),:);
Currents_middle{i} = Current_middle;
Treated_middle{i} = Current_middle(Cut+1:end) %- mean(Current_middle(Cut+1:end));
dumb_Fourier_middle{i} = dumbFT(Treated_middle{i},dt,dw);
FRFFT_middle{i} = FrFFT(Treated_middle{i},dt,dw);
end

%% Plot Middle Currents

figure1 = figure('Name',['Middle Currents, Nc = ',num2str(Ncs),' Nb = ',num2str(Nbs)],'Color',[1 1 1]);

for i = 1:length(U);
    plot(Time,Treated_middle{i},'DisplayName',['U = ',num2str(U(i))],'color',col(i,:))
    %plot(TotalTime,Currents_middle{i},'DisplayName',['U = ',num2str(U(i))],'color',col(i,:))
    hold on
end
hold off
legend('show')
xlabel('Time');
ylabel('Midsection Current');

%% Plot Middle Fouriers

figure1 = figure('Name','Middle Currents, Nc = 50 Nb = 75','Color',[1 1 1]);

for i = 1:length(U);
    %plot(Omega,real(dumb_Fourier_middle{i}),'DisplayName',['DMBFT U = ',num2str(U(i))],'color',col(i,:))
    hold on
    plot(Frequencies,abs(FRFFT_middle{i}),'DisplayName',['FRFFT U = ',num2str(U(i))],'color',col(i,:))
end
hold off
legend('show')
xlabel('Frequencies');
ylabel('Midsection Fourier');
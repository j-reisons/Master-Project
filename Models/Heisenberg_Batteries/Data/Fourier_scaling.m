clc
close all
clear all
%%
J= 1.2;
U = 1.2;

N = [20,50];
Ncs = N;
Nbs = [150,150];%1.5*Ncs;

Cut = 4*Ncs;
Steps = Nbs*10-Cut;

dt = 0.05;
dw = 12./Steps;

Totaltime = cell(1,1);
Time = cell(1,1);
Omega = cell(1,1);
Frequencies = cell(1,1);
windows = cell(1,1);
for i = 1:length(N)
    Totaltime{i} = 0:dt:10*Nbs(i)*dt;
    Time{i} = 0:dt:Steps(i)*dt;
    Omega{i} = 0:dw(i):Steps(i)*dw(i);
    Frequencies{i} = Omega{i}/(2*pi);
    n = [0:Steps(i)];
    l_n = length(n);
    %Blackman-Nuttall
    a = [0.3635819,0.4891775,0.1365995,0.0106411];
    windows{i} = a(1) - a(2)*cos(2*pi*n/l_n) + a(3)*cos(4*pi*n/l_n) - a(4)*cos(6*pi*n/l_n);
    %1 - ((n - ((l_n-1)/2))/((l_n-1)/2)).^2; % Welch
    %0.5*(1-cos(2*pi*n/l_n)); %Hann 
    %
    % 
end

D_max = 100;
freesweeps = 2;
cut = 0.8;

tag = '';
%% Filenames et al.

col=hsv(2*length(N));
Fouriers = cell(1,1);
Currents_middle = cell(1,1);
Treated_middle = cell(1,1);
Windowed_middle = cell(1,1);
FRFFT_treated = cell(1,1);
FRFFT_windowed = cell(1,1);

for i = 1:length(N);
    
filename = ['Batteries','_Nc',strrep(num2str(Ncs(i)),'.',','),'_Nb',...
    strrep(num2str(Nbs(i)),'.',','),'_Ub',strrep(num2str(U),'.',',')...
    ,'_Uc',strrep(num2str(U),'.',','),'_','dt',strrep(num2str(dt),'.',',')...
    ,'_Dmax',num2str(D_max),'_f',num2str(freesweeps),'_',tag,'.mat'];

load(filename);

Current_middle = Currents(round(Nbs(i) + (Ncs(i)/2)),:);
Currents_middle{i} = Current_middle;
Treated_middle{i} = Current_middle(Cut(i)+1:end);% - mean(Current_middle(Cut(i)+1:end));
Windowed_middle{i} = Treated_middle{i}.*windows{i};

FRFFT_treated{i} = FrFFT(Treated_middle{i},dt,dw(i));
FRFFT_windowed{i} = FrFFT(Windowed_middle{i},dt,dw(i));

end

%% Plot Middle Currents

figure1 = figure('Name',['Middle Currents, Nc = ',num2str(Ncs),' Nb = ',num2str(Nbs)],'Color',[1 1 1]);

for i = 1:length(N);
    plot(Time{i},Treated_middle{i},'DisplayName',['N = ',num2str(N(i))],'color',col(i,:))
    %plot(Totaltime{i},Currents_middle{i},'DisplayName',['U = ',num2str(U(i))],'color',col(i,:))
    hold on
end
hold off
legend('show')
xlabel('Time');
ylabel('Midsection Current');

%% Plot Middle Fouriers

figure1 = figure('Name','Middle Currents, Nc = 50 Nb = 75','Color',[1 1 1]);

for i = 1:length(N);
    %plot(Omega,real(dumb_Fourier_middle{i}),'DisplayName',['DMBFT U = ',num2str(U(i))],'color',col(i,:))
    hold on
    semilogy(Frequencies{i},abs(FRFFT_treated{i}),'DisplayName',['treated N = ',num2str(N(i))],'color',col(2*(i-1)+1,:))
    semilogy(Frequencies{i},abs(FRFFT_windowed{i}),'DisplayName',['windowed N = ',num2str(N(i))],'color',col(2*i,:))
end
hold off
legend('show')
xlabel('Frequencies');
ylabel('Midsection Fourier');
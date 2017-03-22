clc
close all
clear all
%%
J= 1;
U_b = 1.1;
U_c = 0;

N = 20:10:50;
dt = 0.05;

D_max = 400;

tag = 'oneU';
%% Filenames et al.

Charges = cell(1);
dC_dts = cell(1);
Interfaces_left = cell(1);
Interfaces_right = cell(1);
Fidelities_cell = cell(1);
Times = cell(1);

i = 1;

col=hsv(2*length(N));

for i = 1:length(N);
    
T = N(i)/(2*J);
    
filename = ['Batteries','_N',strrep(num2str(N(i)),'.',',') ,'_Ub',strrep(num2str(U_b),'.',',')...
    ,'_Uc',strrep(num2str(U_c),'.',','),'_','T',strrep(num2str(T),'.',','),'_','Dmax',num2str(D_max)...
    ,'_',tag,'.mat'];

% filename = ['Batteries','_N',strrep(num2str(N(i)),'.',',') ,'_Ub',strrep(num2str(U_b),'.',',')...
%     ,'_Uc',strrep(num2str(U_c),'.',','),'_','dt',strrep(num2str(dt),'.',','),'_','Dmax',num2str(D_max)...
%     ,'_',tag,'.mat'];

% filename = ['BatteriesXstate','_N',strrep(num2str(N(i)),'.',',') ,'_Ub',strrep(num2str(U_b),'.',',')...
%     ,'_Uc',strrep(num2str(U_c),'.',','),'_','dt',strrep(num2str(dt),'.',','),'_','Dmax',num2str(D_max)...
%     ,'_',tag,'.mat'];

load(filename);

% Magnetization charge
Charge = sum(Magnetizations(1:N(i),:))/N(i);
dC_dt = -N(i)*[0,diff(Charge)]./dt;
Charges{i} = Charge;
dC_dts{i} = dC_dt;

% Interface current left
Interface_left = Currents(N(i),:);
Interfaces_left{i} = Interface_left;

% Interface current right
Interface_right = Currents(2*N(i),:);
Interfaces_right{i} = Interface_right;

% Fidelities
Fidelities_cell{i} = [0,Fidelities];

% Rescaled times
Time = linspace(0,0.5,(10*N(i))+1);
Times{i} = Time;
end

%% Charge plot
figure1 = figure('Name',['Battery charges for Ub = ',num2str(U_b),' Uc = ',num2str(U_c)],'Color',[1 1 1]);

for i = 1:length(N)
    plot(Times{i},Charges{i},'DisplayName',['N = ',num2str(N(i))],'color',col(2*i,:))
    hold on
end
hold off

legend('show')
xlabel('Jt / N');
ylabel('Charge');

%% Current plot
figure2 = figure('Name',['Interface current for Ub = ',num2str(U_b),' Uc = ',num2str(U_c)],'Color',[1 1 1]);

for i = 1:length(N)
    plot(Times{i},Interfaces_left{i},'DisplayName',['Left N = ',num2str(N(i))],'color',col(2*(i-1) + 1,:))
    hold on
    %plot(Times{i},Interfaces_right{i},'DisplayName',['Right N = ',num2str(N(i))],'color',col(i,:))
    plot(Times{i},dC_dts{i},'DisplayName',['dCdT N = ',num2str(N(i))],'color',col(2*i,:))
end
hold off

legend('show')
xlabel('Jt / N');
ylabel('Interface Current');

%% Fidelity plot

figure3 = figure('Name',['Fidelities for Ub = ',num2str(U_b),' Uc = ',num2str(U_c), ' Dmax = ',num2str(D_max)],'Color',[1 1 1]);

for i = 1:length(N)
    semilogy(Times{i},Fidelities_cell{i},'DisplayName',['N = ',num2str(N(i))],'color',col(2*(i-1) + 1,:))
    hold on
end
hold off

legend('show')
xlabel('Jt / N');
ylabel('Fidelities');

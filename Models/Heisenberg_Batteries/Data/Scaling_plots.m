 clc
 %close all
 clear all
%%
J= 1;
U_b = 0;
U_c = 0.5;

N = 20:10:50;
Ncs = N;
Nbs = 1.5*Ncs;
dt = 0.05;

D_max = 100;
freesweeps = 2;

tag = '';
%% Filenames et al.

Charges = cell(1);
dC_dts = cell(1);
Interfaces_left = cell(1);
Interfaces_right = cell(1);
Currents_middle = cell(1);
Fidelities_cell = cell(1);
Times = cell(1);

col=hsv(2*length(N));

for i = 1:length(N);
    
T = N(i)/(2*J);
    
% filename = ['Batteries','_N',strrep(num2str(N(i)),'.',',') ,'_Ub',strrep(num2str(U_b),'.',',')...
%     ,'_Uc',strrep(num2str(U_c),'.',','),'_','T',strrep(num2str(T),'.',','),'_','Dmax',num2str(D_max)...
%     ,'_',tag,'.mat'];

% filename = ['Batteries','_N',strrep(num2str(N(i)),'.',',') ,'_Ub',strrep(num2str(U_b),'.',',')...
%     ,'_Uc',strrep(num2str(U_c),'.',','),'_','dt',strrep(num2str(dt),'.',','),'_','Dmax',num2str(D_max)...
%     ,'_',tag,'.mat'];

% filename = ['BatteriesXstate','_N',strrep(num2str(N(i)),'.',',') ,'_Ub',strrep(num2str(U_b),'.',',')...
%     ,'_Uc',strrep(num2str(U_c),'.',','),'_','dt',strrep(num2str(dt),'.',','),'_','Dmax',num2str(D_max)...
%     ,'_',tag,'.mat'];

filename = ['Batteries','_Nc',strrep(num2str(Ncs(i)),'.',','),'_Nb',...
    strrep(num2str(Nbs(i)),'.',','),'_Ub',strrep(num2str(U_b),'.',',')...
    ,'_Uc',strrep(num2str(U_c),'.',','),'_','dt',strrep(num2str(dt),'.',',')...
    ,'_Dmax',num2str(D_max),'_f',num2str(freesweeps),'_',tag,'.mat'];

load(filename);

% Magnetization charge
Charge = sum(Magnetizations(1:Nbs(i),:))/Nbs(i);
dC_dt = -Nbs(i)*[0,diff(Charge)]./dt;
Charges{i} = Charge;
dC_dts{i} = dC_dt;

% Interface current left
Interface_left = Currents(Nbs(i),:);
Interfaces_left{i} = Interface_left;

% Interface current left
Current_middle = Currents(round(Nbs(i) + (Ncs(i)/2)),:);
Currents_middle{i} = Current_middle;

% Interface current right
Interface_right = Currents(Ncs(i) + Nbs(i),:);
Interfaces_right{i} = Interface_right;

% Fidelities
Fidelities_cell{i} = [0,Fidelities];

% Rescaled times
Time = linspace(0,0.5,(10*Nbs(i))+1);
Times{i} = Time;
end

%% Charge plot
% figure1 = figure('Name',['Battery charges for Ub = ',num2str(U_b),' Uc = ',num2str(U_c)],'Color',[1 1 1]);
% 
% for i = 1:length(N)
%     plot(Times{i},Charges{i},'DisplayName',['N = ',num2str(N(i))],'color',col(2*i,:))
%     hold on
% end
% hold off
% 
% legend('show')
% xlabel('Jt / Nb');
% ylabel('Charge');

%% Current plot
figure2 = figure('Name',['Interface current for Ub = ',num2str(U_b),' Uc = ',num2str(U_c)],'Color',[1 1 1]);

for i = 1:length(N)
    plot(Times{i},Interfaces_left{i},'DisplayName',['Left N = ',num2str(N(i))],'color',col(2*(i-1) + 1,:))
    hold on
    %plot(Times{i},Currents_middle{i},'DisplayName',['Right N = ',num2str(N(i))],'color',col(i,:))
    %plot(Times{i},Interfaces_right{i},'DisplayName',['Right N = ',num2str(N(i))],'color',col(i,:))
    %plot(Times{i},dC_dts{i},'DisplayName',['dCdT N = ',num2str(N(i))],'color',col(2*i,:))
end
hold off

legend('show')
xlabel('Jt / Nb');
ylabel('Interface Current');

%% Fidelity plot

% figure3 = figure('Name',['Fidelities for Ub = ',num2str(U_b),' Uc = ',num2str(U_c), ' Dmax = ',num2str(D_max)],'Color',[1 1 1]);
% 
% for i = 1:length(N)
%     semilogy(Times{i},Fidelities_cell{i},'DisplayName',['N = ',num2str(N(i))],'color',col(2*(i-1) + 1,:))
%     hold on
% end
% hold off
% 
% legend('show')
% xlabel('Jt / Nb');
% ylabel('Fidelities');


% %% Usual plots
% figure
% imagesc(Currents)
% %%
% figure
% imagesc(Magnetizations)
% %%
% %figure
% % plot(Magnetizations(30:50,end))
% % hold on
% %%
% s = size(Magnetizations);
% dM_dT = diff(Magnetizations,1,2);
% dM_dT = [zeros(s(1),1),dM_dT];
% Currents_dM_dT = Currents;
% Currents_dM_dT(1,:) = - dM_dT(1,:);
% for i = 2:s(1)-1
%     Currents_dM_dT(i,:) = -dM_dT(i,:) + Currents_dM_dT(i-1,:);
% end
% Currents_dM_dT = Currents_dM_dT/dt;
% diff = (Currents_dM_dT - Currents);
% diffn = diff;
% diffn(:) = abs(diff(:));
% diffn = diffn(:,50:s(2));
% %%
% figure
% imagesc(diffn)



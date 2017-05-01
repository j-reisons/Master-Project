 clc
 %close all
 clear all
%%
J= 1;
U_b = [1.2879:0.0707:1.7121];
U_c = [1.7121:-0.0707:1.2879];
offdiag = [-0.3:0.1:0.3];

N = 50;
Ncs = N;
Nbs = 1.5*Ncs;
dt = 0.05;

D_max = 500;
freesweeps = 2;
tag = '';

Time = linspace(0,0.5,(10*Nbs)+1);
col=hsv(length(offdiag));

cuts = [0.1,0.45] ; % Interval on which we fit the data
cut_indices = floor(cuts*Nbs/(J*dt));
Cut_Time = Time(cut_indices(1):cut_indices(2)) - Time(cut_indices(1)) + Time(2);

%%
Interfaces_left = cell(1);
Cuts_left = cell(1);

for i = 1:length(offdiag)
    filename = ['Batteries','_Nc',strrep(num2str(Ncs),'.',','),'_Nb',...
        strrep(num2str(Nbs),'.',','),'_Ub',strrep(num2str(U_b(i)),'.',',')...
        ,'_Uc',strrep(num2str(U_c(i)),'.',','),'_','dt',strrep(num2str(dt),'.',',')...
        ,'_Dmax',num2str(D_max),'_f',num2str(freesweeps),'_',tag,'.mat'];
    
    load(filename);
    
    % Interface current left
    Interface_left = Currents(Nbs,:);
    Interface_left(:) = max(Interface_left(:),eps); % Making it non-negative
    Interfaces_left{i} = Interface_left;
    
    % Data which we actually fit
    Cuts_left{i} = Interface_left(cut_indices(1):cut_indices(2));
end

%% Entire time plot
figure1 = figure('Name','Interface current for N = 50, offdiag = -0.3:0.1:0.3','Color',[1 1 1]);
for i = 1 : length(offdiag)
plot(Time,Interfaces_left{i},'DisplayName',['Offdiagonal = ',num2str(offdiag(i))],'color',col(i,:))
hold on
end
hold off

legend('show')
xlabel('Jt / Nb');
ylabel('Interface Current');

%%
figure2 = figure('Name','Fitted data','Color',[1 1 1]);
for i = 1 : length(offdiag)
plot(Cut_Time,Cuts_left{i},'DisplayName',['Offdiagonal = ',num2str(offdiag(i))],'color',col(i,:))
hold on
end
hold off

legend('show')
xlabel('Jt / Nb');
ylabel('Interface Current');

%% Decay rate Fit results

Fitresults  = [-43.82,-37.78,-24.87,-15.95,-18.26,-17.23,-15.32];
Fitresults = -Fitresults;

figure_fitresults = figure('Name','Decay rates','Color',[1 1 1]);
plot(offdiag,Fitresults);
xlabel('Offdiagonal');
ylabel('Decay rate');

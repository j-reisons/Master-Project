function run_battery_Xstate(Nc,Nb,U_c,U_b,dt,D_max,freesweeps,tag)

J = 1;
d = 2;

D_start = round(D_max/4);

Nc
Nb
U_c
U_b
D_max
dt
alpha = 0.8
freesweeps

T = Nb/(2*J);
steps = round(T/dt)

comp_error = 1E-7;

filename = ['BatteriesX','_Nc',strrep(num2str(Nc),'.',','),'_Nb',strrep(num2str(Nb),'.',',')...
    ,'_Ub',strrep(num2str(U_b),'.',',')...
    ,'_Uc',strrep(num2str(U_c),'.',','),'_dt',strrep(num2str(dt),'.',','),'_Dmax',num2str(D_max)...
    ,'_f',num2str(freesweeps),'_',tag,'.mat'];

%% Initial state preparation
State = cell(1,2*Nb + Nc);

Up = zeros(1,1,2);
Up(1,1,1) = 1;
Down = zeros(1,1,2);
Down(1,1,2) = 1;
Xsite = zeros(1,1,2);

Xsite(1,1,1) = 1;
Xsite(1,1,2) = 1;


for i=1:Nb
    State{i} = Up;
    State{Nb + Nc + i} = Down;
end
for i=1:Nc
    State{Nb + i} = Xsite;
end



%%
[U_even_dt,U_odd_dt] = Heisenberg_Batteries_U(Nc,Nb,J,U_b,U_c,dt);
[U_even_half,U_odd_half] = Heisenberg_Batteries_U(Nc,Nb,J,U_b,U_c,dt/2);
U = compressMPO(U_odd_half,U_even_dt,U_odd_half);

Converging_accuracies = cell(1,steps);
Magnetizations = zeros(Nc + 2*Nb,steps+1);
Currents = zeros(Nc + 2*Nb - 1,steps+1);
Fidelities = zeros(1,steps);
sweeps = zeros(1,steps);

S_Z = [1,0;0,-1];
S_X =[0,1;1,0];
S_Y = [0,-1i;1i,0];

S_Z_mpo = cell(1,1);
S_Z_1 = reshape(S_Z,[1,1,d,d]);
S_Z_mpo{1} = S_Z_1;

Q_mpo = cell(1,2);  
Q_1 = zeros(1,2,d,d);
Q_1(1,1,:,:) = S_X*((2*J)^0.5);
Q_1(1,2,:,:) = S_Y*((2*J)^0.5);
Q_mpo{1} = Q_1;
Q_2 = zeros(2,1,d,d);
Q_2(1,1,:,:) = S_Y*((2*J)^0.5);
Q_2(2,1,:,:) = - S_X*((2*J)^0.5);
Q_mpo{2} = Q_2;

%%
State = sweep(State,1);
canon = 1;

for i = 1:steps
    i
    evaluations = Canon_evaluator(State,canon,S_Z_mpo,Q_mpo);
    Magnetizations(:,i) = real(evaluations{1});
    Currents(:,i) = real(evaluations{2});
    
    [State,canon,acc,sw] = cheap_apply_compress(State,U,comp_error,D_max,alpha,freesweeps);
    
    Fidelities(i) = acc(end);
    Converging_accuracies{i} = acc;
    sweeps(i) = sw;
end

evaluations = Canon_evaluator(State,canon,S_Z_mpo,Q_mpo);
Magnetizations(:,steps+1) = real(evaluations{1});
Currents(:,steps+1) = real(evaluations{2});


save(filename,'Magnetizations','Currents','Fidelities','Converging_accuracies','sweeps')
end


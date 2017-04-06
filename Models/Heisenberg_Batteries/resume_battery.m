function resume_battery(Nc,Nb,U_c,U_b,dt,D_max,freesweeps,tag)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

canon = 0;
alpha = 0;
saveperiod = 20;

save_file = make_filename(Nc,Nb,U_b,U_c,dt,D_max,freesweeps,[tag,'_save']);

load(save_file);

for i = i_save:1:steps
    i
    evaluations = Canon_evaluator(State,canon,S_Z_mpo,Q_mpo);
    Magnetizations(:,i) = real(evaluations{1});
    Currents(:,i) = real(evaluations{2});
    
    if toggle
        [State,canon,acc,sw] = cheap_apply_compress_noSVD(State,U,comp_error,alpha,freesweeps);
    else
        [State,canon,acc,sw] = cheap_apply_compress(State,U,comp_error,D_max,alpha,freesweeps);
        if acc > comp_error
            State = Inflate_mps(State,D_max);
            State = sweep(State,canon);
            toggle = 1;
        end
    end
    
    Fidelities(i) = acc(end);
    Converging_accuracies{i} = acc;
    sweeps(i) = sw;
    
    if mod(i,saveperiod) == 0
        i_save = i;
        save(save_file,'Magnetizations','Currents','Fidelities',...
            'Converging_accuracies','sweeps','Nc','Nb','i_save','State',...
            'steps','toggle','U','alpha','S_Z_mpo','Q_mpo','canon','comp_error',...
            'saveperiod');
    end
        
end

evaluations = Canon_evaluator(State,canon,S_Z_mpo,Q_mpo);
Magnetizations(:,steps+1) = real(evaluations{1});
Currents(:,steps+1) = real(evaluations{2});

filename = make_filename(Nc,Nb,U_b,U_c,dt,D_max,freesweeps,[tag,'_resumed']);
save(filename,'Magnetizations','Currents','Fidelities','Converging_accuracies','sweeps','Nc','Nb')

end

function filename = make_filename(Nc,Nb,U_b,U_c,dt,D_max,freesweeps,tag)
filename = ['Batteries','_Nc',strrep(num2str(Nc),'.',','),'_Nb',strrep(num2str(Nb),'.',',')...
    ,'_Ub',strrep(num2str(U_b),'.',',')...
    ,'_Uc',strrep(num2str(U_c),'.',','),'_dt',strrep(num2str(dt),'.',','),'_Dmax',num2str(D_max)...
    ,'_f',num2str(freesweeps),'_',tag,'.mat'];
end

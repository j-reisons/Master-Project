function [mps_out,canon,accs,sweeps] = cheap_apply_compress_noSVD(mps_in,mpo,tol,alpha,freesweeps)
% Memory efficient combination of apply and Iter_comp

%%
N = length(mps_in);

mps_out = sweep(mps_in,-1);

canon = -1;

%%

criterion = 0;
sweeps = 0;

acc = 1;
accs = [];
accs(sweeps + 1) = acc;

R = cell(1,N);
R{N} = 1;

L = cell(1,N);
L{1} = 1;

% Generating all R's
for j = N-1 : -1 : 1
    M_tilde_dag = permute(conj(mps_out{j+1}),[2,1,3]);
    R{j} = contract(R{j+1},2,M_tilde_dag,1);
    
    O_S = apply_site(mpo,mps_in,j+1);
    R{j} = contract(O_S,[2,3],R{j},[1,3]);
    %R{j} = contract(mps_in{j+1},[2,3],R{j},[1,3]);
end

while criterion < 1 && acc > tol
    
    % L -> R sweep
    for i = 1:N
        
        O_S = apply_site(mpo,mps_in,i);
        L_OS = contract(L{i},2,O_S,1);
        %work = contract(L{i},2,mps_in{i},1);
        L_R_OS = contract(L_OS,2,R{i},1);
        mps_out{i} = permute(L_R_OS,[1,3,2]);
        mps_out = L_can(mps_out,i);
        
        %Update L
        M_tilde_dag = permute(conj(mps_out{i}),[2,1,3]);
        L{i+1} = contract(M_tilde_dag,[2,3],L_OS,[1,3]);
        
    end
    sweeps = sweeps+1;
    
    % Stopping criteria
    acc = norm(1 - real(L{N+1}));
    accs(sweeps+1) = acc;
    criterion = (acc/accs(2))*(alpha^(freesweeps - sweeps));
    if criterion > 1 || acc < tol
        canon = 1;
        break
    end
    
    % R -> L sweep
    
    for i = N : -1 : 2
        
        O_S = apply_site(mpo,mps_in,i);
        L_OS = contract(L{i},2,O_S,1);
        %work = contract(L{i},2,mps_in{i},1);
        L_R_OS = contract(L_OS,2,R{i},1);
        mps_out{i} = permute(L_R_OS,[1,3,2]);
        mps_out = R_can(mps_out,i);
        
        %Update R
        M_tilde_dag = permute(conj(mps_out{i}),[2,1,3]);
        R{i-1} = contract(R{i},2,M_tilde_dag,1);
        R{i-1} = contract(O_S,[2,3],R{i-1},[1,3]);
        %R{i-1} = contract(mps_in{i},[2,3],R{i-1},[1,3]);

    end
    sweeps = sweeps+1;
    
    mps_out = R_can(mps_out,1);
    err = contract(R{1},2,conj(mps_out{1}),2);
    O_S = apply_site(mpo,mps_in,1);
    err = contract(err,[1,3],O_S,[2,3]);
    %err = contract(err,[1,3],mps_in{1},[2,3]);

    % Stopping criteria
    acc = norm(1 - real(err));
    accs(sweeps+1) = acc;
    criterion = (acc/accs(2))*(alpha^(freesweeps - sweeps));
    canon = -1;
end
end

function result = apply_site(mpo,mps,site)
result = contract(mps{site},3,mpo{site},4);
R_size = size(result);
result = reshape(permute(result,[3 1 4 2 5]),[R_size(1)*R_size(3),R_size(2)*R_size(4),R_size(5)]);
end


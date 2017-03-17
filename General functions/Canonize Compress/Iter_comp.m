function [mps_out,canon,accs,sweeps] = Iter_comp(mps_in,tol,D_max,alpha,freesweeps)
% Iteratively compresses MPS to D_max and tol. alpha and freesweeps are
% convergence criteria
%
% Stops if tolerance is reached, or if iterations no longer provide
% enough increase.
%
%
% This compression algorithm scales better with bond dimension than
% SVD compression (canonize - compress)
% Outlined in Schollwock 4.5.2 p.46

s = size(mps_in{1});
N = length(mps_in);

mps_out = sweep(mps_in,-1,tol,D_max);
canon = -1;

criterion = 0;
sweeps = 0;

acc = norm(1 - real(braket(mps_in,mps_out)));
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
    R{j} = contract(mps_in{j+1},[2,3],R{j},[1,3]);
end

while criterion < 1 && acc > tol
    
    % L -> R sweep
    for i = 1:N
        
        work = contract(L{i},2,mps_in{i},1);
        work = contract(work,2,R{i},1);
        mps_out{i} = permute(work,[1,3,2]);
        mps_out = L_can(mps_out,i);
        
        %Update L
        M_tilde_dag = permute(conj(mps_out{i}),[2,1,3]);
        L{i+1} = contract(L{i},2,mps_in{i},1);
        L{i+1} = contract(M_tilde_dag,[2,3],L{i+1},[1,3]);
        
    end
    sweeps = sweeps+1;
    
    % Stopping criteria
    prev = acc;
    acc = norm(1 - real(L{N+1}));
    accs(sweeps+1) = acc;
    criterion = (acc/accs(2))*(alpha^(freesweeps - sweeps));
    if criterion > 1 || acc < tol
        canon = 1;
        break
    end
    
    % R -> L sweep
    
    for i = N : -1 : 2
        
        work = contract(L{i},2,mps_in{i},1);
        work = contract(work,2,R{i},1);
        mps_out{i} = permute(work,[1,3,2]);
        mps_out = R_can(mps_out,i);
        
        %Update R
        M_tilde_dag = permute(conj(mps_out{i}),[2,1,3]);
        R{i-1} = contract(R{i},2,M_tilde_dag,1);
        R{i-1} = contract(mps_in{i},[2,3],R{i-1},[1,3]);
    end
    sweeps = sweeps+1;
    
    mps_out = R_can(mps_out,1);
    err = contract(R{1},2,conj(mps_out{1}),2);
    err = contract(err,[1,3],mps_in{1},[2,3]);
    
    % Stopping criteria
    prev = acc;
    acc = norm(1 - real(err));
    accs(sweeps+1) = acc;
    criterion = (acc/accs(2))*(alpha^(freesweeps - sweeps));
    canon = -1;
end
end


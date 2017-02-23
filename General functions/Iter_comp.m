function [mps_out] = Iter_comp(mps_in,D_max,tolerance)
% Compresses MPS to bond dimension D_max iteratively to specified tolerance
% Stops if tolerance is reached, or if iterations no longer provide
% increase in precision

% Function ended up not being used in the code because it did not improve 
% significantly over SVD compression results. Because of that, and
% since it's quite long and that I feel quite lazy this Saturday morning I
% won't comment it.
% Outlined in Schollwock 4.5.2 p.46

s = size(mps_in{1});
N = length(mps_in);

mps_out = sweep(mps_in,1);
mps_out = sweep(mps_out,-1,D_max);

err = zeros(1,5);
acc = 1;
decider = -1;
sweeps = 0;

while decider < -0.1 && acc > tolerance 
	sweeps = sweeps +1;
    R = cell(1,N);
    
    % Generating all R's
    R{N} = 1;
    for j = N-1 : -1 : 1
        M_tilde_dag = permute(conj(mps_out{j+1}),[2,1,3]);
        R{j} = contract(R{j+1},2,M_tilde_dag,1);
        R{j} = contract(mps_in{j+1},[2,3],R{j},[1,3]);
    end
    
    L = 1;
    % L -> R sweep
    for i = 1:N
        
        work = contract(L,2,mps_in{i},1);
        work = contract(work,2,R{i},1);
        mps_out{i} = permute(work,[1,3,2]);
        mps_out = L_can(mps_out,i);
        
        %Update L
        M_tilde_dag = permute(conj(mps_out{i}),[2,1,3]);
        L = contract(L,2,mps_in{i},1);
        L = contract(M_tilde_dag,[2,3],L,[1,3]);
    end
    
    % R -> L sweep
    
    % Generating all L's
    L = cell(1,N);
    L{1} = 1;
    for j = 2 : N
        M_tilde_dag = permute(conj(mps_out{j-1}),[2,1,3]);
        L{j} = contract(L{j-1},2,mps_in{j-1},1);
        L{j} = contract(M_tilde_dag,[2,3],L{j},[1,3]);
    end
    
    R = 1;
    
    for i = N:-1:1
        
        work = contract(L{i},2,mps_in{i},1);
        work = contract(work,2,R,1);
        mps_out{i} = permute(work,[1,3,2]);
        mps_out = R_can(mps_out,i);
        
        %Update R
        M_tilde_dag = permute(conj(mps_out{i}),[2,1,3]);
        R = contract(R,2,M_tilde_dag,1);
        R = contract(mps_in{i},[2,3],R,[1,3]);
    end
    acc = 1 - R;
    
    if sweeps < 5
        err(sweeps) = acc;
    else
        mvavg = mean(err);
        decider = (acc - mvavg)/mvavg;
        mvavg(1) = [];
        mvavg(5) = acc;
    end
    
    
end

end


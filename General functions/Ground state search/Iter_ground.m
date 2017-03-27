function [ground_MPS,energy] = Iter_ground(MPO,MPS,tol)
% Iteratively computes ground state of MPO using MPS as initial guess. Initial guess MUST BE RIGHT CANONIZED (-1)

N = length(MPS);

R = cell(1,N);
R{N} = 1;
L = cell(1,N);
L{1} = 1;
energy = inf;

%Generating all R-expressions

for j = N-1 : -1 : 1
    R{j} = contract(R{j+1},1,MPS{j+1},2);
    if (j == N-1)
        s = size(R{j});
        R{j} = reshape(R{j},[1,1,s(2),s(3)]);
    end
    R{j} = contract(R{j},[1,4],MPO{j+1},[2,4]);
    R{j} = contract(R{j},[1,4],conj(MPS{j+1}),[2,3]);
end

not_eigenstateness = 1;
criterion = 0;

%% L ---> R sweep
while not_eigenstateness > tol && criterion < 0.8
    
    for j = 1:N-1
        s = size(MPS{j});
        opts.v0 = reshape(MPS{j},[s(1)*s(2)*s(3),1]);
        fun = HV(R{j},L{j},MPO{j});
        [M,energy] = eigs(fun,s(1)*s(2)*s(3),1,'SR',opts);
        
        M = reshape(M,[s(1),s(2),s(3)]);
        MPS{j} = M;
        MPS = L_can(MPS,j);
        
        L{j+1} = contract(L{j},1,MPS{j},1);
        if j == 1
            s = size(L{j+1});
            L{j+1} = reshape(L{j+1},[1,s(1),s(2),s(3)]);
        end
        L{j+1} = contract(L{j+1},[1,4],MPO{j},[1,4]);
        L{j+1} = contract(L{j+1},[1,4],conj(MPS{j}),[1,3]);
    end
    
    %% L <--- R sweep
    for j = N:-1:2
        s = size(MPS{j});
        opts.v0 = reshape(MPS{j},[s(1)*s(2)*s(3),1]);
        fun = HV(R{j},L{j},MPO{j});
        [M,energy] = eigs(fun,s(1)*s(2)*s(3),1,'SR',opts);
        
        M = reshape(M,[s(1),s(2),s(3)]);
        MPS{j} = M;
        MPS = R_can(MPS,j);
        
        R{j-1} = contract(R{j},1,MPS{j},2);
        if (j == N)
            s = size(R{j-1});
            R{j-1} = reshape(R{j-1},[1,1,s(2),s(3)]);
        end
        R{j-1} = contract(R{j-1},[1,4],MPO{j},[2,4]);
        R{j-1} = contract(R{j-1},[1,4],conj(MPS{j}),[2,3]);
    end
    
    
    %% Stopping criterions : good enough or not converging fast enough
    HMPS = apply(MPO,MPS);
    HH = braket(HMPS,HMPS);
    H_H = braket(MPS,HMPS)^2;
    previous = not_eigenstateness;
    not_eigenstateness = (HH - H_H)/H_H;
    criterion = not_eigenstateness / previous;
    
end

MPS{1} = R_can(MPS,1);
ground_MPS = MPS;
end

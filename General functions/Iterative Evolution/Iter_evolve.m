function [MPS] = Iter_evolve(MPS,Hamiltonian,dt)
%ITER_EVOLVE INPUT MPS MUST BE RIGHT CANONIZED (-1)
%Based on Unifying time evolution and optimization with mps paper. 

dt = dt/2; % We perform 2 sweeps, evolving by dt/2 at each for a second order integrator

MPO = Hamiltonian;
N = length(MPS);

Right = cell(1,N);
Right{N} = 1;
Left = cell(1,N);
Left{1} = 1;

%Generating all Right-expressions

for j = N-1 : -1 : 1
    Right{j} = contract(Right{j+1},1,MPS{j+1},2);
    if (j == N-1)
        s = size(Right{j});
        Right{j} = reshape(Right{j},[1,1,s(2),s(3)]);
    end
    Right{j} = contract(Right{j},[1,4],MPO{j+1},[2,4]);
    Right{j} = contract(Right{j},[1,4],conj(MPS{j+1}),[2,3]);
end


%% L ---> R sweep

for j = 1:N-1
    
    s = size(MPS{j});
    M0 = reshape(MPS{j},[s(1)*s(2)*s(3),1]);
    
    % We don't want to write or exponentiate a d^2 * D^4 matrix, so we just
    % implement it as a function and Runge-Kutta it up
    
    fun = HV(Right{j},Left{j},MPO{j});
    k1 = -1i*fun(M0);
    k2 = -1i*fun(M0 + (dt/2)*k1);
    k3 = -1i*fun(M0 + (dt/2)*k2);
    k4 = -1i*fun(M0 + dt*k3);
    M = M0 + (dt/6)*(k1 + 2*k2 + 2*k3 + k4);
    
    M = reshape(M,[s(1),s(2),s(3)]);  
    % Here I copypaste the L_can function, since I never intented to need
    % the residue and I do not want to have additional complexity in it
    % just to run this algorithm
    
    work = M;
    s_w = size(work);
    work = permute(work,[1 3 2]);
    work = reshape(work,[s_w(3)*s_w(1),s_w(2)]);
    [Q,R] = qr(work,0);
    U = Q;
    SV = R;
    s_u = size(U);
    A = reshape(U,[s_w(1) s_w(3) s_u(2)]);
    A = permute(A,[1 3 2]);
    MPS{j} = A;
    
    % Now we must update the Left terms and backwards evolve the residue SV
    
    Left{j+1} = contract(Left{j},1,MPS{j},1);
    if j == 1
        s = size(Left{j+1});
        Left{j+1} = reshape(Left{j+1},[1,s(1),s(2),s(3)]);
    end
    Left{j+1} = contract(Left{j+1},[1,4],MPO{j},[1,4]);
    Left{j+1} = contract(Left{j+1},[1,4],conj(MPS{j}),[1,3]);
    
    res0 = SV;
    fun = KV(Right{j},Left{j+1});
    k1 = 1i*fun(res0);
    k2 = 1i*fun(res0 + (dt/2)*k1);
    k3 = 1i*fun(res0 + (dt/2)*k2);
    k4 = 1i*fun(res0 + dt*k3);
    res = res0 + (dt/6)*(k1 + 2*k2 + 2*k3 + k4);
    
    % Back to L_Can copy paste code
    if j ~= N
        MPS{j+1} = contract(res,2,MPS{j+1},1);
    else
        MPS{j} = MPS{j}*sign(SV);
    end
end


%% L <--- R sweep

for j = N:-1:2
    s = size(MPS{j});
    M0 = reshape(MPS{j},[s(1)*s(2)*s(3),1]);
    
    % We don't want to write or exponentiate a d^2 * D^4 matrix, so we just
    % implement it as a function and Runge-Kutta it up
    
    fun = HV(Right{j},Left{j},MPO{j});
    k1 = -1i*fun(M0);
    k2 = -1i*fun(M0 + (dt/2)*k1);
    k3 = -1i*fun(M0 + (dt/2)*k2);
    k4 = -1i*fun(M0 + dt*k3);
    M = M0 + (dt/6)*(k1 + 2*k2 + 2*k3 + k4);
    
    M = reshape(M,[s(1),s(2),s(3)]);

    MPS{j} = M;
    % Here I copypaste the R_can function, since I never intented to need
    % the residue and I do not want to have additional complexity in it
    % just to run this algorithm
    
    work = M;
    s_w = size(work);
    work = reshape(work,[s_w(1),s_w(2)*s_w(3)]);
    [Q,R] = qr(work',0);
    US = R';
    V_dag = Q';
    s_v = size(V_dag);
    B = reshape(V_dag,[s_v(1),s_w(2),s_w(3)]);
    MPS{j} = B;
    
    % Now we must update the Right terms and backwards evolve the residue US
    
    Right{j-1} = contract(Right{j},1,MPS{j},2);
    if (j == N)
        s = size(Right{j-1});
        Right{j-1} = reshape(Right{j-1},[1,1,s(2),s(3)]);
    end
    Right{j-1} = contract(Right{j-1},[1,4],MPO{j},[2,4]);
    Right{j-1} = contract(Right{j-1},[1,4],conj(MPS{j}),[2,3]);
    
    res0 = US;
    fun = KV(Right{j-1},Left{j});
    k1 = 1i*fun(res0);
    k2 = 1i*fun(res0 + (dt/2)*k1);
    k3 = 1i*fun(res0 + (dt/2)*k2);
    k4 = 1i*fun(res0 + dt*k3);
    res = res0 + (dt/6)*(k1 + 2*k2 + 2*k3 + k4);
    
    % R_can copy paste
    if j ~=1
        MPS{j-1} = contract(MPS{j-1},2,res,1);
        MPS{j-1} = permute(MPS{j-1},[1,3,2]);
    else
        MPS{j} = MPS{j}*sign(res);
    end
    
end

end


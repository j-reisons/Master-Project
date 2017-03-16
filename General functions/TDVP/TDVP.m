function [MPS] = TDVP(MPS,Hamiltonian,dt)
%ITER_EVOLVE INPUT MPS MUST BE RIGHT CANONIZED (-1)
%Based on Unifying time evolution and optimization with mps paper. 

half_dt = dt/2; % We perform 2 sweeps, evolving by dt/2 at each for a second order integrator

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
    
    M0 = MPS{j};
    fun = minus_i_HM(Right{j},Left{j},MPO{j});
    M = RK4_step(M0,fun,half_dt);
    
    % Here I copypaste the L_can function, since I never intented to need
    % the residue and I do not want to have additional complexity in it
    % just to run this algorithm
    
    [A,res0] = L_can_res(M);
    MPS{j} = A;
    
    % Now we must update the Left terms and backwards evolve the residue SV
    
    Left{j+1} = contract(Left{j},1,MPS{j},1);
    if j == 1
        s = size(Left{j+1});
        Left{j+1} = reshape(Left{j+1},[1,s(1),s(2),s(3)]);
    end
    Left{j+1} = contract(Left{j+1},[1,4],MPO{j},[1,4]);
    Left{j+1} = contract(Left{j+1},[1,4],conj(MPS{j}),[1,3]);
    
    
    fun = i_KV(Right{j},Left{j+1});
    res = RK4_step(res0,fun,half_dt);
    
    % Back to L_Can copy paste code
    if j ~= N
        MPS{j+1} = contract(res,2,MPS{j+1},1);
    else
        MPS{j} = MPS{j}*sign(res);
    end
end


%% L <--- R sweep

for j = N:-1:2
    
    M0 = MPS{j};
    fun = minus_i_HM(Right{j},Left{j},MPO{j});
    M = RK4_step(M0,fun,half_dt);
    
    [B,res0] = R_can_res(M);
    MPS{j} = B;
    
    % Now we must update the Right terms and backwards evolve the residue US
    
    Right{j-1} = contract(Right{j},1,MPS{j},2);
    if (j == N)
        s = size(Right{j-1});
        Right{j-1} = reshape(Right{j-1},[1,1,s(2),s(3)]);
    end
    Right{j-1} = contract(Right{j-1},[1,4],MPO{j},[2,4]);
    Right{j-1} = contract(Right{j-1},[1,4],conj(MPS{j}),[2,3]);
    
    fun = i_KV(Right{j-1},Left{j});
    res = RK4_step(res0,fun,half_dt);
    
    % R_can copy paste
    if j ~=1
        MPS{j-1} = contract(MPS{j-1},2,res,1);
        MPS{j-1} = permute(MPS{j-1},[1,3,2]);
    else
        MPS{j} = MPS{j}*sign(res);
    end
    
end
end


    function V = RK4_step(V0,fun,dt)
        k1 = fun(V0);
        k2 = fun(V0 + (dt/2)*k1);
        k3 = fun(V0 + (dt/2)*k2);
        k4 = fun(V0 + dt*k3);
        V = V0 + (dt/6)*(k1 + 2*k2 + 2*k3 + k4);
    end

    function [A,res] = L_can_res(M)
        work = M;
        s_w = size(work);
        work = permute(work,[1 3 2]);
        work = reshape(work,[s_w(3)*s_w(1),s_w(2)]);
        [Q,R] = qr(work,0);
        U = Q;
        SV = R;
        res = SV;
        s_u = size(U);
        A = reshape(U,[s_w(1) s_w(3) s_u(2)]);
        A = permute(A,[1 3 2]);
    end

    
    function [B,res] = R_can_res(M)
        work = M;
        s_w = size(work);
        work = reshape(work,[s_w(1),s_w(2)*s_w(3)]);
        [Q,R] = qr(work',0);
        US = R';
        res = US;
        V_dag = Q';
        s_v = size(V_dag);
        B = reshape(V_dag,[s_v(1),s_w(2),s_w(3)]);
    end
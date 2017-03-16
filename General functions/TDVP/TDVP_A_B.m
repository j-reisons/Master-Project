function [MPS] = TDVP_A_B(MPS,Hamiltonian,dt)
%Refactoring TDVP function a bit, fix first and last site evolution

half_dt = dt/2; % We perform 2 sweeps, evolving by dt/2 at each for a second order integrator

MPO = Hamiltonian;
N = length(MPS);

%Generating all Right and Left expressions

[Left,Right] = Gen_LR(MPS,MPO);

%% L ---> R sweep
for j = 1:N
    
    M0 = MPS{j};
    fun = minus_i_HM(Right{j},Left{j},MPO{j});
    M = RK4_step(M0,fun,half_dt);
    
    [A,res] = L_can_res(M);
    MPS{j} = A;
    
    if j ~= N
        Left_evo = contract(Left{j},1,MPS{j},1);
        if j == 1
            s = size(Left_evo);
            Left_evo = reshape(Left_evo,[1,s(1),s(2),s(3)]);
        end
        Left_evo = contract(Left_evo,[1,4],MPO{j},[1,4]);
        Left_evo = contract(Left_evo,[1,4],conj(MPS{j}),[1,3]);
        
        fun = i_KV(Right{j},Left_evo);
        res = RK4_step(res,fun,half_dt);
        
        % Throwing the residue to the right
        MPS{j+1} = contract(res,2,MPS{j+1},1);
        
    else
        MPS{j} = MPS{j}*sign(res);
    end
end

[Left,Right] = Gen_LR(MPS,MPO);
%% L <--- R sweep

for j = N:-1:1
    
    M0 = MPS{j};
    fun = minus_i_HM(Right{j},Left{j},MPO{j});
    M = RK4_step(M0,fun,half_dt);
    
    [B,res] = R_can_res(M);
    MPS{j} = B;
    
    if j ~= 1
        Right_evo = contract(Right{j},1,MPS{j},2);
        if (j == N)
            s = size(Right_evo);
            Right_evo = reshape(Right_evo,[1,1,s(2),s(3)]);
        end  
        Right_evo = contract(Right_evo,[1,4],MPO{j},[2,4]);
        Right_evo = contract(Right_evo,[1,4],conj(MPS{j}),[2,3]);
        
        fun = i_KV(Right_evo,Left{j});
        res = RK4_step(res,fun,half_dt);
        
        % Throwing the residue to the left
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
    
    function [Left,Right] = Gen_LR(MPS,MPO)
    N = length(MPS);
    Right = cell(1,N);
    Right{N} = 1;
    Left = cell(1,N);
    Left{1} = 1;

    MPS = sweep(MPS,-1);
    
    for j = N-1 : -1 : 1
        Right{j} = contract(Right{j+1},1,MPS{j+1},2);
        if (j == N-1)
            s = size(Right{j});
            Right{j} = reshape(Right{j},[1,1,s(2),s(3)]);
        end
        Right{j} = contract(Right{j},[1,4],MPO{j+1},[2,4]);
        Right{j} = contract(Right{j},[1,4],conj(MPS{j+1}),[2,3]);
    end
    
    MPS = sweep(MPS,1);
    
    for j = 1 : N - 1
        Left{j+1} = contract(Left{j},1,MPS{j},1);
        if j == 1
            s = size(Left{j+1});
            Left{j+1} = reshape(Left{j+1},[1,s(1),s(2),s(3)]);
        end
        Left{j+1} = contract(Left{j+1},[1,4],MPO{j},[1,4]);
        Left{j+1} = contract(Left{j+1},[1,4],conj(MPS{j}),[1,3]);
    end  
    
    end
function [ground_state] = Heisenberg_ground_itime(N,U,J,D)
%HEISENBERG_GROUND_ITIME Summary of this function goes here
%   Detailed explanation goes here

MPS = random_mps(N,D,2);

dt = -0.1*1i;

[U_even_dt,U_odd_dt] = Heisenberg_U(N,J,U,dt);
[U_even_half,U_odd_half] = Heisenberg_U(N,J,U,dt/2);
Ham = Heisenberg_H(N,J,U);

U = compressMPO(U_odd_half,U_even_dt,U_odd_half);

energy = braket(MPS,apply(Ham,MPS));
for i = 1:10
    MPS = apply(U,MPS);
    tic
    MPS = Iter_comp(MPS,D,1e-10);
    toc
end
prev_energy = energy;
energy = braket(MPS,apply(Ham,MPS));
criterion = norm((prev_energy - energy)/energy);

while criterion > 5E-5
    
for i = 1:10
    MPS = apply(U,MPS);
    MPS = Iter_comp(MPS,D,1e-10);
end

prev_energy = energy;
energy = braket(MPS,apply(Ham,MPS))
criterion = norm((prev_energy - energy)/energy);
end

ground_state = MPS;
end


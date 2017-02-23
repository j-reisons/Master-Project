clc
close all
clear all
%% Parameters
N = 20;
J = 1;
U = 0;
T = 10;
D_max = 10;
dt = 5*10^-2;
%% State , evolution operators , S_Z operators
State = cell(1,N);

for i = 1:N
    State{i} = reshape([0,1],[1,1,2]);
end

State{1} = reshape([1,0],[1,1,2]);

[U_odd,U_even] = Heisenberg_U(N,J,U,dt);

S_z = [1,0;0,-1];

%%
time = linspace(0,T,floor(T/dt)+1);
Profiles = zeros(length(time),N);

for i = 1 : length(time);
    for j = 1:N
        Sz_State = State;
        Sz_State{j} = contract(State{j},3,S_z,2);
        Profiles(i,j) = real(braket(Sz_State,State));
    end
    State = apply(U_odd,State);
    State = sweep(State,1);
    State = sweep(State,-1,D_max);
    State =  apply(U_even,State);
    State = sweep(State,-1);
    State = sweep(State,1,D_max);
end

%%
%save('N = 9, D = 10','Profiles','N')
%% Animation part
%load('N = 9, D = 10');
%load('N = 9,Comparison');
%%
x = 1:N;
for i = 1:length(Profiles)
y_MPS = Profiles(i,:);
plot(x,y_MPS)
legend('MPS');
axis([0.5 N+0.5 -1.2 1.2])
pause(1/60)
end
    
    
clc
close all
clear all
%% Load something
% load('C:\Users\Juris\Documents\MATLAB\Master Project\Models\Heisenberg_Batteries\Batteries_N10_Ub0_Uc1,1_T5_Dmax100_threeU.mat')
s = size(Magnetizations);
N = s(1)/3;
%% Magnetization charge
Charge = sum(Magnetizations(1:N,:))/N;
Transferred_charge = 1 - Charge;
figure
plot(Charge)
%% Interface current
Interface = Currents(N,:);
figure
plot(Interface);
%%
figure
imagesc(Currents)
%%
figure
imagesc(Magnetizations)
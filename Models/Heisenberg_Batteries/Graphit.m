clc
close all
clear all
%%
load('N=10_J=2_Uc=1_Ub=0_T=5_Dmax=100.mat')

oneU = Magnetizations;

load('Batteries_N10_Ub0_Uc1_T5_Dmax100_testing.mat')

threeU = Magnetizations;

diff = oneU-threeU;
%%
imagesc(oneU)
figure
imagesc(threeU)
figure
imagesc(diff)
%%
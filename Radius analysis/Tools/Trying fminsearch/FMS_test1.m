close all
clear all
clc


R = 80;
T = 20;
N = 100000;


xy = MC_Sim_3DShell(N, R, T);
r = sqrt(xy(:,1).^2 + xy(:,2).^2);
r_hist = 0:1:ceil(R+T/2+1);

figure;
hist(r,r_hist);
n = hist(r,r_hist);

h_Chi2_fun = @(rt0) CalcChi2(rt0,r_hist,n);




%%
r_init = 90;
t_init = 30;

rt_init = [r_init t_init];
tic
[rt_opt,fval] = fminsearch(h_Chi2_fun,rt_init);
toc

[rt_opt,fval]
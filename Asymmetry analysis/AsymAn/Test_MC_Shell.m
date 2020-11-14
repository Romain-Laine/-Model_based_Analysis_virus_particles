clear all
close all
clc

N_sim = 20;
r_hist = linspace(0,25,10);

r0 = 20
n_sample = 100

figure;

for i = 1:10
xy_sim = MC_Sim_3DShell(N_sim,r0,0);
r_sim = sqrt(xy_sim(:,1).^2 + xy_sim(:,2).^2);
n_Shell = histc(r_sim,r_hist);
% n_Shell = n_sample*n_Shell/sum(n_Shell);

hold on
plot(r_hist,n_Shell)
end


xy_sim = MC_Sim_3DShell(N_sim,r0,0);

figure;
plot(xy_sim(:,1),xy_sim(:,2),'+')
axis equal

Theta = atan2(xy_sim(:,2),xy_sim(:,1))/pi;

Theta_hist = -1:0.1:1;
hTheta = histc(Theta,Theta_hist);

n_av = N_sim/(size(Theta_hist,2)-1);

figure;
plot(Theta_hist,hTheta,Theta_hist,ones(1,size(Theta_hist,2))*n_av,'r')
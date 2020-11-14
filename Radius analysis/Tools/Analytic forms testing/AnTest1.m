close all
clear all
clc

r = 60;
t = 10;

ErrorLoc = 20;


% N = 100;
% rho = linspace(0,r+t,N);

rho_Step = 0.1;
Max_rho = r + t/2 + 2*ErrorLoc + 1;
rho = 0:rho_Step:Max_rho ;

V_1 = sqrt((r+t/2)^2 - rho.^2).*rho;
V_1(rho > r+t/2) = 0;

V_2 = sqrt((r-t/2)^2 - rho.^2).*rho;
V_2(rho > r-t/2) = 0;

V_shell = V_1 - V_2;

figure;
plot(rho,V_1,rho,V_2, rho, V_shell)

V_shell = V_shell/sum(V_shell);

%%

x = (-4*ErrorLoc):rho_Step:(4*ErrorLoc);
G_ErrorLoc = exp(-x.^2/(2*ErrorLoc^2));
G_ErrorLoc = G_ErrorLoc/sum(G_ErrorLoc);

figure;
plot(x,G_ErrorLoc);

V_Error = conv(cat(2,zeros(1,(size(x,1)-1)/2),V_shell),G_ErrorLoc,'same');

V_1e = sqrt((r+t/2)^2 - rho.^2);
V_1e(rho > r+t/2) = 0;

V_2e = sqrt((r-t/2)^2 - rho.^2);
V_2e(rho > r-t/2) = 0;

V_shell_E = V_1 - V_2;
V_Error2 = conv(V_shell_E,G_ErrorLoc,'same');

V_Error2 = V_Error2.*rho;

figure;
plot(rho,V_Error2)


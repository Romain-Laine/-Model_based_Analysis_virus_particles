close all
clear all
clc


tic
r = 78;
t = 30;

Linker_size = 10;
Linker_width = 2;

LocError = 20;
LocError2D = LocError*sqrt(2);

Sim_Step = 1;   % nm


% ----------------------------------------------------------------
x_prot = -(r+t/2):Sim_Step:(r+t/2);
[X_prot,Y_prot] = meshgrid(x_prot);
Z_prot = real(sqrt((r+t/2)^2 - X_prot.^2 - Y_prot.^2)) - real(sqrt((r-t/2)^2 - X_prot.^2 - Y_prot.^2));

figure;
mesh(X_prot,Y_prot,Z_prot);
%%
% ---------------
x_linker = -(Linker_size + Linker_width/2):Sim_Step:(Linker_size + Linker_width/2);
[X_linker,Y_linker] = meshgrid(x_linker);
Z_linker = real(sqrt((Linker_size+Linker_width/2)^2 - X_linker.^2 - Y_linker.^2)) - real(sqrt((Linker_size-Linker_width/2/2)^2 - X_linker.^2 - Y_linker.^2));

figure;
mesh(X_linker,Y_linker,Z_linker);


% ---------------
x_LocError = -(4*LocError):Sim_Step:(4*LocError);
[X_LocError,Y_LocError] = meshgrid(x_LocError);
Z_LocError = exp(-(X_LocError.^2 + Y_LocError.^2)/(2*LocError2D^2));

figure;
mesh(X_LocError,Y_LocError,Z_LocError);

% Convolution
Z_tot = conv2(Z_LocError,conv2(Z_prot,Z_linker));
toc

figure;
imshow(Z_tot,[])

%%
n = size(Z_tot,1);
MyRow = Z_tot((n-1)/2+1,(n-1)/2+1:end);
% MyRow = Z_tot((n-1)/2+1,:);
r = 0:Sim_Step:Sim_Step*(n-1)/2;
rHisto = abs(r).*MyRow;

figure;
plot(r,rHisto )
xlim([0, (n-1)/2*Sim_Step])






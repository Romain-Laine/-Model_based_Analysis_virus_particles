%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code simulates a single-color STORM image (localization) for clustering, and uses Ripley's K analysis
% Dr Romain Laine, rfl30@cam.ac.uk, Laser Analytics Group 2015-11-03
%
% - Simulate density of points with clustering (Gaussian clustering)
% - Display a STORM image
% - Running K clustering algorithm
% - Calculting CI bound using Monte Carlo simulation of CSR
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all
clc

Fov = 20000; % nm  10 um is roughly 64*64 pixels
SR_image_pixelSize = 32; % nm

% K-ripley analysis
r_step = 10;
Analysis_window = 500; % nm (no point going further than 1um)
% CI calculation
N_repeat_CI = 100;
p = 0.01; % if p = 0.01 it is equal to 99% CI, it's probably good enough

% --- Cluster parameters ---
N_clusters = 500;
N_loc_per_cluster = 1;
Cluster_size = 100;  % the cluster size is defined as 2*standard deviatrion of the Gaussian distribution

Cluster_density = N_clusters /(Fov/1000)^2;
disp(['Red cluster density: ', num2str(Cluster_density),' clusters/um^2']);


% ------------------------------------------------------------------------------------
% Generate cluster dataset

[x_c,y_c] = Generate_random_XY( Fov, N_clusters );
[ X_cluster, Y_cluster] = XYcentre2Cluster( x_c, y_c, N_loc_per_cluster, Cluster_size, Fov);

%%
SR_imageSize = ceil(Fov/SR_image_pixelSize)+2;
[ SR_image_rS ] = SRdisplay_rapidSTORM( cat(2,X_cluster,Y_cluster), SR_image_pixelSize, SR_imageSize );
ColorImage_rS = Grey2Color(SR_image_rS/max(SR_image_rS(:)));

figure('Color','white','name','Two-Color SR image');
imshow(ColorImage_rS);
title([num2str(SR_image_pixelSize),' nm / pixel'])

%%

[ r_hist, ~, K, ~, D ] = CalcRipleyK( X_cluster, Y_cluster, X_cluster, Y_cluster, Fov, Fov^2, Analysis_window, r_step );
H = sqrt(K/pi) - r_hist;  % H(r) = L(r) - r
[ ~, ~, ~, CI_upperbound, CI_lowerbound ] = RipleyK_CIcalc( Fov, D, D, Analysis_window, r_step, N_repeat_CI, p );
% The CI is calculated on a FOV that is the avergae of all of the FOV

%%
figure('Color','white','name','L(r)-r Ripley''s function');
plot(r_hist,H)
hold on
plot(r_hist, CI_upperbound, '--r')
hold on
plot(r_hist, CI_lowerbound, '--r')
xlim([0 Analysis_window])
xlabel('r (nm)')
ylabel('Ripley''s L(r)-r function')
grid on




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code simulates 2-color STORM images for clustering, and uses Ripley's K analysis
% Dr Romain Laine, rfl30@cam.ac.uk, Laser Analytics Group 2015-11-03
%
% - Simulate density of points with clustering as well as 2 color co-clustering
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
r_step = 20;
Analysis_window = 500; % nm (no point going further than 1um)
% CI calculation
N_repeat_CI = 10000; % the minimum is 100 !
p = 0.01; % if p = 0.01 it is equal to 99% CI, it's probably good enough

% --- Cluster 1 ---
N_clusters_1 = 50;
N_loc_per_cluster_1 = 5;
Cluster_size_1 = 50;

% --- Cluster 2 ---
N_clusters_2 = 25;
N_loc_per_cluster_2 = 5;
Cluster_size_2 = 30;

% --- Dual-color parameters ---
Correlated = 1;
f = 50;   % in %, fraction of loc. having dual color
Loc_distance = 200;
Loc_distance_flexibility = 1000; % +/- 0.5 x Loc_distance_flexibility

% ------------------------------------------------------------------------------------
% Generate cluster dataset 1

[x_c1,y_c1] = Generate_random_XY( Fov, N_clusters_1 );
[ X_cluster1, Y_cluster1] = XYcentre2Cluster( x_c1, y_c1, N_loc_per_cluster_1, Cluster_size_1, Fov);

% Generate cluster dataset 2
if Correlated == 0
    [x_c2, y_c2] = Generate_random_XY( Fov, N_clusters_2 );
elseif Correlated == 1
    [x_c2, y_c2] = GenerateCorrelated_XY( x_c1, y_c1, Loc_distance, Loc_distance_flexibility, f);
end
[ X_cluster2, Y_cluster2] = XYcentre2Cluster( x_c2, y_c2, N_loc_per_cluster_2, Cluster_size_2, Fov);

RedCluster_density = length(X_cluster1)/(Fov/1000)^2;
GreenCluster_density = length(X_cluster2)/(Fov/1000)^2;

disp(['Red cluster density: ', num2str(RedCluster_density),' clusters/um^2']);
disp(['Green cluster density: ', num2str(GreenCluster_density),' clusters/um^2']);


%%
SR_imageSize = ceil(Fov/SR_image_pixelSize)+2;
[ SR_image_rS1 ] = SRdisplay_rapidSTORM( cat(2,X_cluster1,Y_cluster1), SR_image_pixelSize, SR_imageSize );
ColorImage_rS1 = Grey2Color(SR_image_rS1/max(SR_image_rS1(:)));

[ SR_image_rS2 ] = SRdisplay_rapidSTORM( cat(2,X_cluster2,Y_cluster2), SR_image_pixelSize, SR_imageSize );
ColorImage_rS2 = Grey2SummerColor( SR_image_rS2/max(SR_image_rS2(:)));

figure('Color','white','name','Two-Color SR image');
imshow(ColorImage_rS1 + ColorImage_rS2);
title([num2str(SR_image_pixelSize),' nm / pixel'])
%%

[ r_hist, ~, K12, ~, D12 ] = CalcRipleyK( X_cluster1, Y_cluster1, X_cluster2, Y_cluster2, Fov, Fov^2, Analysis_window, r_step );
H12 = sqrt(K12/pi) - r_hist;  % H(r) = L(r) - r
[ ~     , ~, K21, ~, D21 ] = CalcRipleyK( X_cluster2, Y_cluster2, X_cluster1, Y_cluster1, Fov, Fov^2, Analysis_window, r_step );
H21 = sqrt(K21/pi) - r_hist;  % H(r) = L(r) - r

%%

[ ~, ~, ~, CI_upperbound12, CI_lowerbound12 ] = RipleyK_CIcalc( Fov, D12, D21, Analysis_window, r_step, N_repeat_CI, p );
[ ~, ~, ~, CI_upperbound21, CI_lowerbound21 ] = RipleyK_CIcalc( Fov, D21, D12, Analysis_window, r_step, N_repeat_CI, p );
% The CI is calculated on a FOV that is the avergae of all of the FOV

%%

figure('Color','white','name','L(r)-r functions', 'Units','normalized','Position',[0.1 0.1 0.75 0.4]);
subplot(1,2,1)
plot(r_hist, H12)
hold on
plot(r_hist, CI_upperbound12, '--r')
hold on
plot(r_hist, CI_lowerbound12, '--r')
xlim([0 Analysis_window])
xlabel('r (nm)')
ylabel('Ripley''s L(r)-r function')
grid on

subplot(1,2,2)
plot(r_hist, H21)
hold on
plot(r_hist, CI_upperbound21, '--r')
hold on
plot(r_hist, CI_lowerbound21, '--r')
xlim([0 Analysis_window])
xlabel('r (nm)')
ylabel('Ripley''s L(r)-r function')
grid on




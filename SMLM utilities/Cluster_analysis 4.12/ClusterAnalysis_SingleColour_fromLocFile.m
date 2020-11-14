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

disp('Single Color cluster analysis 4.10 from Localization files')
% Remove the warning message saying the image is too big for display
WarningID = 'images:initSize:adjustingMag';
warning('off',WarningID);

CameraPixelSize = 157.2; % nm
Number_of_pixels = 256;

Fov = CameraPixelSize*Number_of_pixels; % nm  10 um is roughly 64*64 pixels
SR_image_pixelSize = 32; % nm  for displaying the SR image

File_token = '_RC';
area_token = 'a';   % file name is built as "area_token"+area number+"File_token"

% K-ripley analysis
r_step = 10;
Analysis_window = 500; % nm (no point going further than 1um)
% CI calculation
N_repeat_CI = 100;
p = 0.01; % if p = 0.01 it is equal to 99% CI, it's probably good enough
CI_calculation = 0;

% Display
n = 8; % number of bits for color images

DefaultPath = 'C:\Users\rfl30\DATA raw\2015_12_15 Janin''s data for clustering\for Romain_clusteranalysis_150128\dopaminergic neurons 17122015\control\';
% [FileName,FolderName] = uigetfile('*.txt','Select a localization file',DefaultPath);
% disp(['Loading ',FolderName,FileName]);

[ PathName, FileNames ] = GetListofLocFile( DefaultPath, File_token, area_token );
N_files = size(FileNames,1);
disp('Folder chosen: ');
disp(PathName);
disp(['Number of files found: ', num2str(N_files)]);

%%

X_all = cell(N_files,2);
Y_all = cell(N_files,2);

Area_all = zeros(N_files,1);
Fov_all = zeros(N_files,2);

for i = 1:N_files
    tic
    disp('------------------------------------------------');
    disp('------------------------------------------------');
    disp(['Opening: ',FileNames{i}]);
    
    LocFile = Read_LocFile( FileNames{i}, 0);
    disp(['Number of localizations: ', num2str(size(LocFile,1))]);
    
    X = LocFile(:,1);
    Y = LocFile(:,2);
    
    SR_imageSize = ceil(Fov/SR_image_pixelSize)+2;
    SR_image_rS = SRdisplay_rapidSTORM( cat(2,X,Y), SR_image_pixelSize, SR_imageSize );
    ColorImage_rS = Grey2Color(SR_image_rS/max(SR_image_rS(:)));
    Fov_all(i,1) = size(SR_image_rS,1)*SR_image_pixelSize; % nm
    Fov_all(i,2) = size(SR_image_rS,2)*SR_image_pixelSize; % nm
    
    figure('Color','white','name',[area_token, num2str(i),File_token,' - Two-Color SR image']);
    imshow(ColorImage_rS);
    title([num2str(SR_image_pixelSize),' nm / pixel'])
    
    Mask = FindAreaOutline_SRimage( SR_image_rS, Analysis_window, SR_image_pixelSize );
    Area = sum(Mask(:))*SR_image_pixelSize^2;   % Area in nm^2
    
    ColorImage(:,:,1) = (2^n-1)*SR_image_rS/max(SR_image_rS(:));
    ColorImage(:,:,3) = (2^n-1)*bwmorph(Mask,'remove');
    figure('Color','white','name',[area_token, num2str(i),File_token,' - Mask']);
    imshow(ColorImage);
    title([num2str(SR_image_pixelSize),' nm / pixel']);
    
    %%
    X_all{i,1} = X;  % this structure is required to accommodate 2 colour clustering
    X_all{i,2} = X;
    
    Y_all{i,1} = Y;
    Y_all{i,2} = Y;
    
    
    Area_all(i) = sum(Mask(:))*SR_image_pixelSize^2;   % Area in nm^2
    
    disp(['Total FOV: ', num2str((Fov_all(i,1)*Fov_all(i,2))/10^6), ' um^2 (',num2str(Fov_all(i,1)/10^3), ' um x ', num2str(Fov_all(i,2)/10^3),' um)']);
    disp(['Analysis area: ', num2str(Area_all(i)/10^6), ' um^2']);
    disp(['Number of localizations: ', num2str(length(X)), ' (local density: ',num2str(10^6*length(X)/Area_all(i)), ' localizations / um^2)']);
    toc
    
end

r_hist = (0:r_step:Analysis_window)';
disp('------------------------------------------------');
disp('------------------------------------------------');
disp('Computing Ripley K function...')
[ ~, H, ~, D ] = CalculateRipleyKH_combineFOV( X_all, Y_all, Area_all, Analysis_window, r_step, Fov_all );


if CI_calculation == 1
    
    disp('------------------------------------------------');
    disp('------------------------------------------------');
    disp('Computing confidence interval...')
    %     FOV_for_CI = [mean(Fov_all(:,1)), mean(Fov_all(:,2))];
    %     FOV_for_CI = [sqrt(mean(Area_all)), sqrt(mean(Area_all))];
    FOV_for_CI = [sqrt(mean(Fov_all(:,1))/mean(Fov_all(:,2))*mean(Area_all)), sqrt(mean(Fov_all(:,2))/mean(Fov_all(:,1))*mean(Area_all))];
    
    [ ~, ~, ~, CI_bound ] = RipleyK_CIcalc_singleColor( FOV_for_CI, D, Analysis_window, r_step, N_repeat_CI, p );
    % The CI is calculated on a FOV that is the avergae of all of the FOV
    
    figure('Color','white','name','L(r)-r Ripley''s function');
    plot(r_hist,H)
    hold on
    plot(r_hist, CI_bound, '--r')
    xlim([0 Analysis_window])
    xlabel('r (nm)')
    ylabel('Ripley''s L(r)-r function')
    grid on
    
else
    figure('Color','white','name','L(r)-r Ripley''s function');
    plot(r_hist,H)
    xlim([0 Analysis_window])
    xlabel('r (nm)')
    ylabel('Ripley''s L(r)-r function')
    grid on
    
end

warning('on',WarningID);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code analyses 2-color STORM images for clustering, by using Ripley's K analysis
% Dr Romain Laine, rfl30@cam.ac.uk, Laser Analytics Group 2015-11-03
%
% - Object segmentation using Otsu thresholding
% - Excluding points for corection of edge effect
% - Running K clustering algorithm
% - Calculting CI bound using Monte Carlo simulation of CSR
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all
clc

disp('Dual Color cluster analysis 4.12 from Images')
% Remove the warning message saying the image is too big for display
WarningID = 'images:initSize:adjustingMag';
warning('off',WarningID);

% ------------------------------------------------------------------------------------
% dSTORM image display
SR_image_pixelSize = 30; % nm

% Centroids detection parameters
RC_Gamma = 0.2;
RC_Otsu_level = 1;
RC_Min_object_size = 10; % in pixels !

GC_Gamma = 0.2;
GC_Otsu_level = 1;
GC_Min_object_size = 10; % in pixels !

% K-ripley analysis
r_step = 5;
Analysis_window = 500; % nm (no point going further than 1um)

% Display
n = 8; % number of bits for color images

% CI calculation
N_repeat_CI = 10000;
p = 0.01; % if p = 0.01 it is equal to 99% CI, it's probably good enough
CI_calculation = 1;

% Range for averaging
r0 = 100; % nm
r1 = 150; % nm

% ------------------------------------------------------------------------------------
DefaultPath = 'C:\Users\rfl30\DATA raw\dSTORM v1 data\Images from Colin for Clustering\';

% RC_token = '_647_stack-STORMdataImage';
% GC_token = '_568_stack-STORMdataImage';
% area_token = 'area';

RC_token = '_RC';
GC_token = '_GC_registered';
area_token = 'a';

% ------------------------------------------------------------------------------------
t_start = tic;
[ PathName, FileNames_RC, FileNames_GC ] = Get2CSRimage( DefaultPath, RC_token, GC_token, area_token );
N_files = size(FileNames_RC,1);
disp('------------------------------------------------');
disp('------------------------------------------------');

disp('Opening files in following folder:');
disp(PathName);
disp(['Number of files found: ', num2str(N_files)]);

% ------------------------------------------------------------------------------------
% Open measured dataset

X_all12 = cell(N_files,2);
Y_all12 = cell(N_files,2);
X_all21 = cell(N_files,2);
Y_all21 = cell(N_files,2);
Area_all = zeros(N_files,1);
Fov_all = zeros(N_files,2);

for i = 1:N_files
    tic
    disp('------------------------------------------------');
    disp('------------------------------------------------');
    
    disp(['Opening: ',FileNames_RC{i}]);
    Image_RC = imread(FileNames_RC{i});
    Fov_all(i,1) = size(Image_RC,2)*SR_image_pixelSize; % nm (here we do not assume that the image is square)
    Fov_all(i,2) = size(Image_RC,1)*SR_image_pixelSize; % nm (here we do not assume that the image is square)
    % X direction is horizontal Y is vertical in MATLAB display
    
    [ X1, Y1, ~] = SRimage2Centroids( Image_RC, SR_image_pixelSize, RC_Gamma, RC_Otsu_level, RC_Min_object_size );
    
    disp('------------------------------------------------');
    disp(['Opening: ',FileNames_GC{i}]);
    Image_GC = imread(FileNames_GC{i});
    [ X2, Y2, ~] = SRimage2Centroids( Image_GC, SR_image_pixelSize, GC_Gamma, GC_Otsu_level, GC_Min_object_size );
    [ X2, Y2 ] = RemoveOutFOV_loc( X2,Y2, Fov_all(i,:) ); % the green color may be 2C registered and therefore we may find localization outside the RC field of view
    %     therefore this step is necessary to get rid of them in the GC only.
    
    SR_image_rS1 = SRdisplay_rapidSTORM( cat(2,X1,Y1), SR_image_pixelSize, max(size(Image_RC)) );   % Red channel
    SR_image_rS2 = SRdisplay_rapidSTORM( cat(2,X2,Y2), SR_image_pixelSize, max(size(Image_GC)) );   % Green channel
    
    Mask = FindAreaOutline( SR_image_rS1 + SR_image_rS2, Analysis_window, SR_image_pixelSize );
    ColorImage = Create2Cimage( SR_image_rS1, SR_image_rS2, n ); % RC and GC in that order
    
    ColorImage(:,:,3) = (2^n-1)*bwmorph(Mask,'remove');
    figure('Color','white','name','Two-Color SR image');
    imshow(ColorImage);
    title([num2str(SR_image_pixelSize),' nm / pixel']);
    
    X_all12{i,1} = X1;
    X_all12{i,2} = X2;
    
    Y_all12{i,1} = Y1;
    Y_all12{i,2} = Y2;
    
    X_all21{i,1} = X2;
    X_all21{i,2} = X1;
    
    Y_all21{i,1} = Y2;
    Y_all21{i,2} = Y1;
    
    Area_all(i) = sum(Mask(:))*SR_image_pixelSize^2;   % Area in nm^2
    
    disp(['Total FOV: ', num2str((Fov_all(i,1)*Fov_all(i,2))/10^6), ' um^2']);
    disp(['Analysis area: ', num2str(Area_all(i)/10^6), ' um^2']);
    disp(['Number of objects in RC: ', num2str(length(X1)), ' (local density: ',num2str(10^6*length(X1)/Area_all(i)), ' centroids / um^2)']);
    disp(['Number of objects in GC: ', num2str(length(X2)), ' (local density: ',num2str(10^6*length(X2)/Area_all(i)), ' centroids / um^2)']);
    toc
    
end


%%

r_hist = (0:r_step:Analysis_window)';

disp('------------------------------------------------');
disp('------------------------------------------------');
disp('Computing Ripley K function...')
[ ~, H12, ~, D12 ] = CalculateRipleyKH_combineFOV( X_all12, Y_all12, Area_all, Analysis_window, r_step, Fov_all );
[ ~, H21, ~, D21 ] = CalculateRipleyKH_combineFOV( X_all21, Y_all21, Area_all, Analysis_window, r_step, Fov_all );

H_av12 = mean(H12((r_hist >= r0) & (r_hist <= r1)));
H_av21 = mean(H21((r_hist >= r0) & (r_hist <= r1)));

disp(['Average H12 (between ', num2str(r0), ' and ', num2str(r1), ' nm): ', num2str(H_av12)]);
disp(['Average H21 (between ', num2str(r0), ' and ', num2str(r1), ' nm): ', num2str(H_av21)]);

if CI_calculation == 1
    disp('------------------------------------------------');
    disp('------------------------------------------------');
    disp('Computing confidence interval...')
    %     FOV_for_CI = [mean(Fov_all(:,1)), mean(Fov_all(:,2))];
    %     FOV_for_CI = [sqrt(mean(Area_all)), sqrt(mean(Area_all))];
    FOV_for_CI = [sqrt(mean(Fov_all(:,1))/mean(Fov_all(:,2))*mean(Area_all)), sqrt(mean(Fov_all(:,2))/mean(Fov_all(:,1))*mean(Area_all))];
    
    [ ~, ~, ~, CI_bound12 ] = RipleyK_CIcalc_dualColor( FOV_for_CI, D12, D21, Analysis_window, r_step, N_repeat_CI, p );
    [ ~, ~, ~, CI_bound21 ] = RipleyK_CIcalc_dualColor( FOV_for_CI, D21, D12, Analysis_window, r_step, N_repeat_CI, p );
    % The CI is calculated on a FOV that is the avergae of all of the FOV
    
    HCI_av12 = mean(CI_bound12((r_hist >= r0) & (r_hist <= r1)));
    HCI_av21 = mean(CI_bound21((r_hist >= r0) & (r_hist <= r1)));
    
    disp(['Average H12 CSR (between ', num2str(r0), ' and ', num2str(r1), ' nm): ', num2str(HCI_av12), '(ratio = ', num2str(H_av12/HCI_av12), ')']);
    disp(['Average H21 CSR (between ', num2str(r0), ' and ', num2str(r1), ' nm): ', num2str(HCI_av21), '(ratio = ', num2str(H_av21/HCI_av21), ')']);
    
end

%% Display the H functions

if CI_calculation == 0
    figure('Color','white','name','L(r)-r functions with CI bound', 'Units','normalized','Position',[0.1 0.1 0.75 0.4]);
    subplot(1,2,1)
    plot(r_hist, H12)
    xlim([0 Analysis_window])
    xlabel('r (nm)')
    ylabel('Ripley''s L(r)-r function')
    grid on
    legend('Data (G on R)', 'Location','Southeast')
    title 'Clustering of the green channel data on the red channel data'
    
    subplot(1,2,2)
    plot(r_hist, H21)
    xlim([0 Analysis_window])
    xlabel('r (nm)')
    ylabel('Ripley''s L(r)-r function')
    grid on
    legend('Data (R on G)', 'Location','Southeast')
    title 'Clustering of the red channel data on the green channel data'
    
    Results = [r_hist H12 H21];
    
elseif CI_calculation == 1
    
    figure('Color','white','name','L(r)-r functions with CI bound', 'Units','normalized','Position',[0.1 0.1 0.75 0.4]);
    subplot(1,2,1)
    plot(r_hist, H12)
    hold on
    plot(r_hist, CI_bound12, '--r')
    xlim([0 Analysis_window])
    xlabel('r (nm)')
    ylabel('Ripley''s L(r)-r function')
    grid on
    legend('Data (G on R)', [num2str(100*(1-p)),'% CI CSR'],'Location','Southeast')
    title 'Clustering of the green channel data on the red channel data'
    
    subplot(1,2,2)
    plot(r_hist, H21)
    hold on
    plot(r_hist, CI_bound21, '--r')
    xlim([0 Analysis_window])
    xlabel('r (nm)')
    ylabel('Ripley''s L(r)-r function')
    grid on
    legend('Data (R on G)', [num2str(100*(1-p)),'% CI CSR'],'Location','Southeast')
    title 'Clustering of the red channel data on the green channel data'
    
    Results = [r_hist H12 CI_bound12 H21 CI_bound21];
end


disp('------------------------------------------------');
toc(t_start);
warning('on',WarningID);




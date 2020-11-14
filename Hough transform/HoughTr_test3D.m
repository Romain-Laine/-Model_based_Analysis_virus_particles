clear all
close all
clc

% DATA = 'Simulated';
DATA = 'Measured';

ProtStructure = 'Shell';
% ProtStructure = 'Ring';

LinkerStructure = 'Shell';
% LinkerStructure = 'Ring';

% Pixel size for density display
PixelSize = 10; % nm

x0 = -60;
y0 = 10;

radVP = 100;            % nm
ThicknessVP = 30;       % nm
Linker_size = 10;       % nm
Linker_thickness = 2;   % nm
LocError = 20;          % nm

% Total number of localizations
N_loc = 500;    % number of localizations (assuming 1 fluorophore / protein)


%------------------------------------------------------------------------------------------------------------------------%
% Generate the protein localization

if strcmp(DATA,'Simulated')
    disp(['Simulated radius and thickness:      ',num2str([radVP,ThicknessVP],4)]);
    disp(['N = ',num2str(N_loc)]);
    tic
    if strcmp(ProtStructure,'Shell')
        Prot_loc = MC_Sim_3DShell(N_loc, radVP, ThicknessVP);
    elseif strcmp(ProtStructure,'Ring')
        Prot_loc = MC_Sim_2DRing(N_loc, radVP, ThicknessVP);
    end
    
    % Add the offset
    Prot_loc(:,1) = Prot_loc(:,1) + x0;
    Prot_loc(:,2) = Prot_loc(:,2) + y0;
    
    % Generate the fluorophore localization
    if strcmp(LinkerStructure,'Shell')
        Fluo_loc = Prot_loc + MC_Sim_3DShell(N_loc, Linker_size, Linker_thickness);
    elseif strcmp(LinkerStructure,'Ring')
        Fluo_loc = Prot_loc + MC_Sim_2DRing(N_loc, Linker_size, Linker_thickness);
    end
    
    % Generate the error on the localization due to precision
    xy_loc_error = MC_Sim_2DLocError(N_loc,LocError);
    
    % Sum it all up to get the localization
    xy = Fluo_loc + xy_loc_error;
    toc
    
    Size = max(abs(x0),abs(y0)) + 2*radVP;
    Size = PixelSize*round(Size/PixelSize);
    
    x_grid = -Size:PixelSize:Size;
    xy_round = PixelSize*round(xy/PixelSize);
    Density_image = zeros(size(x_grid,2),size(x_grid,2));
    
    %%
    for i = 1:size(xy_round,1)
        Density_image((x_grid == xy_round(i,1)),(x_grid == xy_round(i,2))) = Density_image((x_grid == xy_round(i,1)),(x_grid == xy_round(i,2))) + 1;
    end
    
    %------------------------------------------------------------------------------------------------------------------------%
    % Display it
    figure('Color','white');
    plot(xy(:,1),xy(:,2),'+')
    hold on
    plot(Prot_loc(:,1),Prot_loc(:,2),'go')
    hold on
    plot(Fluo_loc(:,1),Fluo_loc(:,2),'rx')
    axis equal
    xlabel 'x (nm)'
    ylabel 'y (nm)'
    grid on
    
    Cmap = colormap(hot(256));
    Density_image = ApplyCmap2gsc(Density_image, Cmap);
    
end


% Using measured dataset
if strcmp(DATA,'Measured')
    % Load dataset
    DefaultPath = 'C:\Users\rfl30\Desktop\Hough transform test images\';
    [FileName,PathName,FilterIndex] = uigetfile([DefaultPath,'*.png']);
    Density_image = imread([PathName, FileName],'png');
end

Density_image_raw = Density_image;


Smoothing = 1;
if Smoothing == 1
    GaussFilter = fspecial('gaussian',[3 3],0.5);
    Density_image = imfilter(Density_image,GaussFilter,'same');
    figure('Color','white');
    subplot(1,2,1)
    imshow(Density_image,[]);
    title 'Raw image'
    subplot(1,2,2)
    imshow(Density_image,[]);
    title 'With Gaussian filter'
elseif Smoothing == 0
    figure('Color','white');
    imshow(Density_image,[]);
    tite 'No smoothing applied'
end


%%
Radius_range = [100 200]; % nm
allowed_radius = round(Radius_range/PixelSize); % in pixels


% allowed_radius = round(radVP/PixelSize*[0.5, 1.5]);
% allowed_radius = round([60 200]/PixelSize);
disp(['Allowed radius range: ',num2str(allowed_radius), ' (in pixels)']);

% Perform the Hough transform
Sensitivity = 0.98;
Method = 'PhaseCode';
% Method = 'TwoStage';
EdgeThresh = 0.2;

[centers, radii] = imfindcircles(Density_image,allowed_radius,'ObjectPolarity', 'bright', 'Sensitivity',Sensitivity,'Method',Method,'EdgeThreshold',EdgeThresh);

n_circles = size(radii,1);
disp(['Number of circles found: ',num2str(n_circles)]);

if strcmp(DATA,'Simulated')
    disp(['Simulated values: ',num2str(radVP),' ',num2str(y0),' ',num2str(x0)]);
end

% Calculate the radius and centre in nm (taking the pixel size into account)
calc_radius = radii*PixelSize;
calc_center = PixelSize*centers - PixelSize*(1+(size(Density_image,2)-1)/2);

disp('Radii (nm):');
disp(num2str(calc_radius));

disp(' ');
disp('Centers (nm):')
disp(num2str(calc_center));

%%


Image_Hough = Density_image;
Raw_Image_Hough = Density_image_raw;

for i = 1:n_circles
    ImCircle = CreateCircleImage(Density_image, centers(i,:), radii(i) );
    
    Raw_Image_Hough(round(centers(i,2)),round(centers(i,1)),2) = 255;
    Image_Hough(round(centers(i,2)),round(centers(i,1)),2) = 255;
    
    Raw_Image_Hough(:,:,3) = Raw_Image_Hough(:,:,3) + uint8(255*ImCircle);
    Image_Hough(:,:,3) = Image_Hough(:,:,3) + uint8(255*ImCircle);
end

figure('Color','white');
subplot(2,2,1)
imshow(Density_image_raw,[]);
title 'Raw image'
subplot(2,2,2)
imshow(Density_image,[]);
title 'With Gaussian filter'
subplot(2,2,3)
imshow(Raw_Image_Hough);
title 'Hough circle detection (raw)'
subplot(2,2,4)
imshow(Image_Hough);
title 'Hough circle detection (smoothed)'



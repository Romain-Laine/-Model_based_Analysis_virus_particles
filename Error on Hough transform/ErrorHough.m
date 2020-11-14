%%%%%%%%%%%%% Estimating the error on the centre of particles using HOugh transform %%%%%%%%%%%%%
%% 2014-06-09
% Romain Laine (rfl30@cam.ac.uk)

%------------------------------------------------------------------------------------------------------------------------%
% Clear out
clear all
close all
clc

%------------------------------------------------------------------------------------------------------------------------%
ProtStructure = 'Shell';
% ProtStructure = 'Ring';

LinkerStructure = 'Shell';
% LinkerStructure = 'Ring';

radVP = 78.9;             % nm
ThicknessVP = 10.6;     % nm
Linker_size = 15;       % nm
Linker_thickness = 5;   % nm
LocError = 15;          % nm

% Total number of localizations
N_loc = 273;    % number of localizations/particle (assuming 1 fluorophore / protein)
PixelSize = 10;
SR_imageSize = 80;   % in pixels  (basically 800 nm / pixelsize)

Display = 0;

%------------------------------------------------------------------------------------------------------------------------%
% Generate the protein localization

disp(['Simulated radius: ',num2str(radVP),' nm - thickness: ',num2str(ThicknessVP),' nm']);
N_sim = 10000;
x0y0 = zeros(N_sim,2);
xHyH = zeros(N_sim,2);
rH = zeros(N_sim,1);

tic
h_wait = waitbar(0,'Please wait...');
for k = 1:N_sim
    waitbar(k/N_sim,h_wait);
    if strcmp(ProtStructure,'Shell')
        Prot_loc = MC_Sim_3DShell(N_loc, radVP, ThicknessVP);
    elseif strcmp(ProtStructure,'Ring')
        Prot_loc = MC_Sim_2DRing(N_loc, radVP, ThicknessVP);
    end
    
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
    
    %------------------------------------------------------------------------------------------------------------------------%
    % Display it without offset
    if Display == 1
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
    end
    
    % Generate the centre position randomly within 200-500 nm central square
    x0 = 300 + rand(1)*200;
    y0 = 300 + rand(1)*200;
    x0y0(k,1) = x0;
    x0y0(k,2) = y0;
    
    xy(:,1) = xy(:,1) + x0;
    xy(:,2) = xy(:,2) + y0;
    
    % Display the intensity map
    if Display == 1
        disp('----------------------------');
        disp('Histogramming...');
        SR_image_hist = SRimage_histogram(xy, SR_imageSize, PixelSize);
        figure('name','Intensity image based on simple histogramming');
        imshow(SR_image_hist,[])
    end
    
    SR_image = SRdisplay_rapidSTORM( xy, PixelSize, SR_imageSize);
    if Display == 1
            disp('rapidSTORM pixelization...');

        figure('name','Intensity image based on rapidSTORM pixelization');
        imshow(SR_image,[])
    end
    
    %% Hough transform
    
    SR_image_raw = SR_image;
    
    % Apply some Gaussian smoothing
    Smoothing = 1;
    if Smoothing == 1
        GaussFilter = fspecial('gaussian',[5 5],0.75);
        SR_image = imfilter(SR_image,GaussFilter,'same');
        
        if Display == 1
                    disp('Gaussian smoothing applied.');

            figure('Color','white','name','Hough transform - Smoothing');
            subplot(1,2,1)
            imshow(SR_image_raw,[]);
            title 'Raw image'
            subplot(1,2,2)
            imshow(SR_image,[]);
            title 'With Gaussian filter'
        end
    elseif Smoothing == 0
        disp('No smoothing applied.');
        if Display == 1
            figure('Color','white','name','Hough transform - Smoothing');
            imshow(SR_image,[]);
            title 'No smoothing applied'
        end
    end
    
    % Good initial values are as following:
    Radius_range = [70 200]; % nm
    Sensitivity = 0.93;
    EdgeThresh = 0.05;
    
    % Convert in pixels
    allowed_radius = round(Radius_range/PixelSize); % in pixels
    
    %% Perform the Hough transform
    
    Method = 'PhaseCode';
    % Method = 'TwoStage';
    if Display == 1
        disp(['Allowed radius range: ',num2str(allowed_radius), ' (in pixels)']);
        disp(['Sensitivity: ',num2str(Sensitivity),' - Edge threshold: ', num2str(EdgeThresh)]);
        disp('Computing Hough transform...');
    end
    [centers, radii, circle_strength] = imfindcircles(SR_image,allowed_radius,'ObjectPolarity', 'bright', 'Sensitivity',Sensitivity,'Method',Method,'EdgeThreshold',EdgeThresh);

    
    % Total number of circles found by the Hough transform
    n_circles = size(radii,1);
    centres_nm = (centers-1)*PixelSize;
    
    if size(centers,1) ~= 1 || size(centers,2) ~= 2
        xHyH(k,:) = [0,0];
        rH(k) = 0;
    else
        xHyH(k,:) = centres_nm;
        rH(k) = radii;
    end
    
    if Display == 1
        [ ImCircle, ImCentres, LogicalIm ] = CreateCircleImage(size(SR_image), centers, radii);
        SR_image_8bit = uint8(255*SR_image/max(SR_image(:)));
        Image_Hough = cat(3,SR_image_8bit,SR_image_8bit,SR_image_8bit);
        Image_Hough(:,:,3) = Image_Hough(:,:,3) + uint8(255*ImCircle);
        Image_Hough(:,:,2) = Image_Hough(:,:,2) + uint8(255*ImCentres);
        
        figure('name','Particle with detected circle');
        imshow(Image_Hough);
    end
end
close(h_wait);
toc

%%
xHyH(rH == 0,:) = [];
x0y0(rH == 0,:) = [];
rH(rH == 0) = [];


%%
dxy = xHyH - x0y0;

[n_dx, dx_hist] = hist(dxy(:,1),50);
[n_dy, dy_hist] = hist(dxy(:,2),50);

[fitresult_x, gof_x] = GaussFit(dx_hist, n_dx);
[fitresult_y, gof_y] = GaussFit(dy_hist, n_dy);

%% Total
[n_dxy, dxy_hist] = hist(cat(1,dxy(:,1),dxy(:,2)),50);
[fitresult_xy, gof_xy] = GaussFit(dxy_hist, n_dxy);

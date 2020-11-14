%%%%%%%%%%%%% Radius analysis %%%%%%%%%%%%%
%% 2014-06-02
% Romain Laine (rfl30@cam.ac.uk)

% This code analyses localization data and fits a virus particle model to it.
% The virus particle model is as follows:
% The protein structure and the structure of the linker can be chosen among several different simulation:
% - Ring (circle with thickness) (2D), circle is a particlar case of ring with thickness = 0
% - Shell (3D), sphere is a particular case of shell with thickness = 0
% NB: the disk and the sphere volume can be obtained by picking thickness = 2x radius

% The analysis is performed by minimization of the residuals between simulated and data radius histogram (projected on 2D, as is the data).
% The minimization is performed in multiple stages:
% - using fminsearch MATLAB function (2x)
% - using non-linear least square fitting using Levenberg-Marquardt algorithm

%------------------------------------------------------------------------------------------------------------------------%
% Clear out
clear all
close all
clc

tStart = tic;
%------------------------------------------------------------------------------------------------------------------------%

% For measured dataset
Default_PathName = 'E:\STORM data\2014_06_03 Radius analysis\';

% When saving the images of the fitting
Save_path = 'C:\Users\rfl30\Desktop\Testing save\';

DATA = 'Simulated';
% DATA = 'Measured';

%------------------------------------------------------------------------------------------------------------------------%
ProtStructure = 'Shell';
% ProtStructure = 'Ring';

LinkerStructure = 'Shell';
% LinkerStructure = 'Ring';

radVP = 80;             % nm
ThicknessVP = 10;       % nm

Linker_size = 10;       % nm
Linker_thickness = 2;   % nm
LocError = 15;          % nm

% Total number of localizations
N_loc = 10000;    % number of localizations (assuming 1 fluorophore / protein)
Display_SRimage = 1;


%------------------------------------------------------------------------------------------------------------------------%
% Pixel size for density display (display only)
PixelSize = 10;

% For rapidSTORM display
ScaleBarSize = 100; % nm
% SR_imageSize = max(xy_rapidS(:))/PixelSize + 5;
SR_imageSize = 40;


%------------------------------------------------------------------------------------------------------------------------%
% Generate the protein localization
disp(['Dataset type: ', DATA]);
if strcmp(DATA,'Simulated')
    disp(['Simulated radius: ',num2str(radVP),' nm - thickness: ',num2str(ThicknessVP),' nm']);
    
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
    % Calculate radius and histograms
    r = sqrt((xy(:,1)).^2 + (xy(:,2)).^2);
    Max_r = radVP + 0.5*ThicknessVP + Linker_size + 2*LocError + 1;
    
    %------------------------------------------------------------------------------------------------------------------------%
    % Display it
    figure('Color','white','name','Simulated localizations');
    plot(xy(:,1),xy(:,2),'+')
    hold on
    plot(Prot_loc(:,1),Prot_loc(:,2),'go')
    hold on
    plot(Fluo_loc(:,1),Fluo_loc(:,2),'rx')
    axis equal
    xlabel 'x (nm)'
    ylabel 'y (nm)'
    grid on
    
    MaxRad = Max_r;
    
end


% Using measured dataset
if strcmp(DATA,'Measured')
    
    % Load dataset
    [ParticleInfo, FolderName] = Multiple_COM(Default_PathName);
    disp(FolderName);
    
    prompt = {'Mininum frame number:','Minimum ADC:','Maximum radius (nm):'};
    dlg_title = 'Input';
    num_lines = 1;
    def = {'0','1750','150'};
    answer = inputdlg(prompt,dlg_title,num_lines,def);
    
    if isempty(answer)
        MinFrame = 0;
        MinADC = 0;
        MaxRad = 40000; % 40 um, it should be enough
    else
        MinFrame = str2double(answer{1});
        MinADC = str2double(answer{2});
        MaxRad = str2double(answer{3});
    end
    
    
    disp('----------------------------');
    disp(['Min. ADC: ',num2str(MinADC),' ADC']);
    disp(['Min. frame: ',num2str(MinFrame),' fr.']);
    disp(['Max. radius: ',num2str(MaxRad),' nm']);
    
    
    xy = ParticleInfo(:,[1,2]);
    LocRad = sqrt((xy(:,1)).^2 + (xy(:,2)).^2);
    Frame = ParticleInfo(:,3);
    ADC = ParticleInfo(:,4);
    ParticleInfo(ADC < MinADC | Frame < MinFrame | LocRad > MaxRad,:) = [];
    
    disp('Fixed parameters:');
    disp(['Linker size: ',num2str(Linker_size),' nm +/- ',num2str(Linker_thickness),' nm']);
    disp(['Localization error: ',num2str(LocError),' nm']);
    
    % Calculate radius and histograms
    xy = ParticleInfo(:,[1,2]);
    r = sqrt((xy(:,1)).^2 + (xy(:,2)).^2);
    Max_r = prctile(r,99.5);
    
    %------------------------------------------------------------------------------------------------------------------------%
    % Display it
    figure('Color','white','name','Localizations image');
    plot(xy(:,1),xy(:,2),'+')
    axis equal
    xlabel 'x (nm)'
    ylabel 'y (nm)'
    grid on
end

%------------------------------------------------------------------------------------------------------------------------%
disp(['Total number of localizations in dataset: ',num2str(size(r,1))]);

% Visualize the data (bin number is not important here)
% r_hist = linspace(0,Max_r,min(round(N_loc/4),round(Max_r)));
r_hist = 0:1:200;


figure('Color','white','name','Radius histogram');
hist(r,r_hist)
xlim([0,Max_r])
xlabel 'Radius (nm)'
ylabel 'Occurences'
title(['Mean radius: ',num2str(mean(r),4),' and STD: ',num2str(std(r),4)])




if Display_SRimage == 1
    % Display the intensity map
    PixelNumber = round(Max_r/PixelSize);
    Size = PixelNumber*PixelSize;
    SRimage_histogram(xy, PixelNumber,Size);
    
    %% rapidSTORM display
    xy_rapidS = xy;
    xy_rapidS(:,1) = xy_rapidS(:,1) + SR_imageSize*PixelSize/2;
    xy_rapidS(:,2) = xy_rapidS(:,2) + SR_imageSize*PixelSize/2;
    
    [ SR_image_rS ] = SRdisplay_rapidSTORM( xy_rapidS, PixelSize, SR_imageSize );
    ColorImage_rS = Grey2Color(SR_image_rS/max(SR_image_rS(:)));
    ColorImage_rS = AddScaleBar( ColorImage_rS, PixelSize, ScaleBarSize);
    figure;
    imshow(ColorImage_rS)
    title(['Scale bar: ',num2str(ScaleBarSize),' nm (',num2str(PixelSize),' nm / pixel)'])
    
    imwrite(ColorImage_rS,[Save_path,'SR_image_r',num2str(radVP),'t',num2str(ThicknessVP),'.png']);
    %
end

disp('----------------------------');

[n_r,r_c] = hist(r,r_hist);
n_r = n_r';
r_c = r_c';


figure('Color','white','name','Radius histogram');
plot(r_c,n_r);



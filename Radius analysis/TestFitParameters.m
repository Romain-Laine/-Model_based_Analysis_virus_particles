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

%------------------------------------------------------------------------------------------------------------------------%

% For measured dataset
Default_PathName = 'E:\STORM data\2014_06_03 Radius analysis\';

% DATA = 'Simulated';
DATA = 'Measured';

%------------------------------------------------------------------------------------------------------------------------%
ProtStructure = 'Shell';
% ProtStructure = 'Ring';

LinkerStructure = 'Shell';
% LinkerStructure = 'Ring';

radVP = 68;             % nm
ThicknessVP = 10;       % nm
Linker_size = 20;       % nm
Linker_thickness = 5;   % nm
LocError = 20;          % nm

% Total number of localizations
N_loc = 25000;    % number of localizations (assuming 1 fluorophore / protein)

%------------------------------------------------------------------------------------------------------------------------%
% Fitting parameters
Weighted = 'on';
% Weighted = 'off';

N_display = 10000000; % number of localization to simulate for the display

% For SR image display only
PixelSize = 10; % nm

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
    
end


% Using measured dataset
if strcmp(DATA,'Measured')
    
    % Load dataset
    [ParticleInfo, FolderName] = Multiple_COM(Default_PathName);
    disp(FolderName);
    
    prompt = {'Mininum frame number:','Minimum ADC:','Maximum radius (nm):'};
    dlg_title = 'Input';
    num_lines = 1;
    def = {'0','0','150'};
    answer = inputdlg(prompt,dlg_title,num_lines,def);
    
    if isempty(answer)
        MinFrame = 0;
        MinADC = 0;
        MinRad = 40000; % 40 um, it should be enough
    else
        MinFrame = str2double(answer{1});
        MinADC = str2double(answer{2});
        MinRad = str2double(answer{3});
    end
    
    
    disp('----------------------------');
    disp(['Min. ADC: ',num2str(MinADC),' ADC']);
    disp(['Min. frame: ',num2str(MinFrame),' fr.']);
    disp(['Max. radius: ',num2str(MinRad),' nm']);
    
    
    xy = ParticleInfo(:,[1,2]);
    LocRad = sqrt((xy(:,1)).^2 + (xy(:,2)).^2);
    Frame = ParticleInfo(:,3);
    ADC = ParticleInfo(:,4);
    ParticleInfo(ADC < MinADC | Frame < MinFrame | LocRad > MinRad,:) = [];

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
disp('----------------------------');


%%
%------------------------------------------------------------------------------------------------------------------------%
% Calculate the data histogram

% Maximum number of bins permissible
n_bin_max = 150;
% Number of bins used for the histogramming
n_bin = min([round(N_loc/4),round(Max_r),n_bin_max]);

% The bins are chosen to be linearly spaced (is there a wiser choice?)
r_hist = linspace(0,Max_r,n_bin);
N_data = hist(r,r_hist);

% Eliminate the last bin as it may contain all the ouliers in measured
% dataset --> skews results
N_data(end) = [];

disp(['Total number of localizations used for histogramming: ',num2str(sum(N_data))]);
disp(['Number of bins in histogram: ',num2str(n_bin)]);



%% Display the results for the optimal parameters
close all

radVP = 66.5;             % nm
ThicknessVP = 38.2;       % nm
Linker_size = 20;       % nm
Linker_thickness = 5;   % nm
LocError = 15;          % nm

disp('----------------------------');
disp('Parameters:');
disp(['Particle radius: ',num2str(radVP),' nm']);
disp(['Particle thickness: ',num2str(Linker_thickness),' nm']);
disp(['Linker size: ',num2str(Linker_size),' nm +/- ',num2str(Linker_thickness),' nm']);
disp(['Localization error: ',num2str(LocError),' nm']);


% initialize sim_param variable
Sim_param = cell(5,1);
Sim_param{1} = ProtStructure;
Sim_param{2} = LinkerStructure;
Sim_param{3} = Linker_size;        % nm
Sim_param{4} = Linker_thickness;   % nm
Sim_param{5} = LocError;           % nm

N_opt = VirusSim(Sim_param, radVP, ThicknessVP, N_display, r_hist);
N_opt = N_opt*sum(N_data)/sum(N_opt);

%%
figure('Color','white','name','Radius histogram - optimal parameters');
plot(r_hist(1:end-1),N_data,'+')
hold on
plot(r_hist(1:end-1),N_opt,'r')
xlabel 'Radius (nm)'
ylabel 'Occurences'
grid on
title(['Radius: ',num2str(radVP,'%3.2f'), ' nm - Thickness: ',num2str(ThicknessVP,'%3.2f'),' nm'])

%%

% Calculate Chi^2

% Weight array
if strcmp(Weighted, 'on')
    W = N_data;
    W(W == 0) = 1;
    W = 1./W;
elseif strcmp(Weighted, 'off')
    W = ones(1,size(N_data,2));
end

% Calculate the degree of freedom for CHi2 calculation
n_freedom = size(N_data,2) - 2 - 1;
Chi2 = (1/n_freedom)*sum(W.*(N_data-N_opt).^2);
disp(['Chi^2: ',num2str(Chi2)]);


%% Display SR image

LocInfo = zeros(size(xy,1),4);
LocInfo(:,1) = xy(:,1) - min(xy(:,1));
LocInfo(:,2) = xy(:,2) - min(xy(:,2));
LocInfo(:,3) = ones(size(xy,1),1);
LocInfo(:,4) = ones(size(xy,1),1);

GreyImage = Loc2SRimage(LocInfo, PixelSize);
ColorImage = Grey2Color(GreyImage);

ScaleBarSize = 100; % nm
ColorImage = AddScaleBar( ColorImage, PixelSize, ScaleBarSize);

figure;
imshow(ColorImage,[])


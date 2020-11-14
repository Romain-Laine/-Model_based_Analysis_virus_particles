%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code analyses the localization files obtained from rapidSTORM to
% estimate the localization precision and resolution
% Dr Romain Laine, rfl30@cam.ac.uk, Laser Analytics Group 2015-11-10
%
% - Choose appropriate camera settings
% - Load localization files
% - Threshold the ones below Min_ADC
% - COmpute the localization precision histogram and statistics
% - The resolution is estimated as the FWHM of the localization
% distribution (e.g. 2*sqrt(2*log(2))*<Loc. precision>)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all
clear all
clc
disp('ADC to localization precision 8.0');

n_file = 7; % Number of files
DefaultPath = 'C:\Users\rfl30\DATA raw\';

AnalyseTracks = 0; % set to 1 if you want to check how long are the emission in frames

% Camera settings;
CamSettings = 5;

if CamSettings == 1 % Settings for the old camera
    Min_ADC = 0; % ADC (as usually set on rapidSTORM)
    PixelSize = 158.4;        % nm
    Lambda = 676;             % nm
    NA = 1.03;
    Bg_noise = 25;            % ADC
    Excess_noise = sqrt(2);   % EMCCD excess noise
    Conversion_factor = 1.76; % ADC/photon
    
elseif CamSettings == 2 % Red channel
    Min_ADC = 850; % ADC (as usually set on rapidSTORM)
    PixelSize = 157.2;        % nm
    Lambda = 676;             % nm
    NA = 1.03;
    Bg_noise = 25;            % ADC
    Excess_noise = sqrt(2);   % EMCCD excess noise
    Conversion_factor = 3.15; % ADC/photon
    
elseif CamSettings == 3 % Green channel
    Min_ADC = 1000; % ADC (as usually set on rapidSTORM)
    PixelSize = 157.2;        % nm
    Lambda = 607;             % nm
    NA = 1.0;
    Bg_noise = 35;            % ADC
    Excess_noise = sqrt(2);   % EMCCD excess noise
    Conversion_factor = 3.15; % ADC/photon
    
elseif CamSettings == 4 % Blue channel STORMv2
    Min_ADC = 6700; % ADC (as usually set on rapidSTORM)
    PixelSize = 118;        % nm
    Lambda = 520;             % nm
    NA = 0.7;
    Bg_noise = 150;            % ADC on high noise camera acquisition from SMTI
    Excess_noise = sqrt(2);   % EMCCD excess noise
    Conversion_factor = 13.6; % ADC/photon
    
elseif CamSettings == 5 % Blue channel STORMv2
    Min_ADC = 10001; % ADC (as usually set on rapidSTORM)
    PixelSize = 117;        % nm
    Lambda = 676;             % nm
    NA = 1;
    Bg_noise = 12;            % ADC on high noise camera acquisition from SMTI
    Excess_noise = sqrt(2);   % EMCCD excess noise
    Conversion_factor = 12.3; % ADC/photon
    
end

disp('----------------');
disp(['Camera settings: ',num2str(CamSettings)]);
disp(['Min ADC: ',num2str(mean(Min_ADC)),' ADC']);
disp(['Background noise: ',num2str(Bg_noise),' ADC']);
disp(['Conversion factor: ',num2str(Conversion_factor),' ADC/photon']);

% calculate the PSF sigma
PSF_sigma = 0.22*Lambda/NA; % assuming Gaussian description of the PSF
Sigma_T = sqrt(PSF_sigma^2 + PixelSize^2/12);
PSF_FWHM = 2*sqrt(2*log(2))*PSF_sigma;
disp(['PSF FWHM: ',num2str(PSF_FWHM), ' nm' ]);

% Compute the Localization precision using Mortensen formulae
LocPrec = @(x) Excess_noise*sqrt(Sigma_T^2*Conversion_factor./x.*(16/9 + 8*pi()*Sigma_T^2*(Bg_noise/Conversion_factor)^2*Conversion_factor./(x*PixelSize^2)));
ADC = 200:10:20000;

figure('color','white','name','Localization precision vs. ADC');
plot(ADC,LocPrec(ADC));
xlabel 'ADC'
ylabel 'Localization precision (nm)'
grid on


%% Load the localization file

if AnalyseTracks == 1
    ParticleInfo = zeros(0,5);
else
    ParticleInfo = zeros(0,4);
end


for i = 1:n_file
    [FileName,FolderName] = uigetfile('*.txt','Select a localization file',DefaultPath);
    disp(['Loading ',FolderName,FileName]);
    DefaultPath = FolderName; % Set the default path to the current folder
    
    LocFile = Read_LocFile( [FolderName, FileName], AnalyseTracks);
    
    
    disp(['Number of localizations: ', num2str(size(LocFile,1))]);
    ParticleInfo = cat(1,ParticleInfo, LocFile );
end

% Remove ADC below the Min ADC
disp('----------------');
ADCs = ParticleInfo(:,4);
Total_number_loc = length(ADCs);
disp(['Total number of localizations for analysis: ', num2str(Total_number_loc)]);
disp(['Removing localization below ', num2str(Min_ADC), ' ADC...']);
ParticleInfo(ADCs <= Min_ADC,:) = [];
ADCs = ParticleInfo(:,4);

Number_loc = length(ADCs);
disp(['Fraction of localization left: ', num2str(Number_loc), ' (', num2str(100*Number_loc/Total_number_loc,'%.1f'), ' %)']);

% Compute the localization precisions from these ADCs
LocPrecision = LocPrec(ADCs);

% Round the precisions to the nearest 0.5 nm
round5LocPrecision = 5*round(round(10*LocPrecision)/5)/10;

% Display results
disp('----------------');
disp(['Number of localizations: ',num2str(length(ADCs))]);
disp('----------------');
disp('Localization precision statistics (Mortensen formulae):')
disp(['Mean: ',num2str(mean(LocPrecision),'%.1f'),' nm']);
disp(['Mode: ',num2str(mode(round5LocPrecision)),' nm']);
disp(['Median: ',num2str(median(round5LocPrecision)),' nm']);
disp(['STD: ',num2str(std(LocPrecision),'%.1f'),' nm']);
disp('----------------');
FWHM_loc_precision = 2*sqrt(2*log(2))*mean(LocPrecision);
disp(['Resolution: ',num2str(FWHM_loc_precision,'%.1f'),' nm (as calculated as FWHM of localization)' ]);

%%

figure('name','Localization precision histogram','Units','normalized', 'Position', [0.15 0.15 0.8 0.5],'color','white');
subplot(1,2,1)
ADCs_hist = 0:100:(max(ADCs)+100);
histogram(ADCs, ADCs_hist);
xlabel 'Counts (ADC)'
ylabel 'Occurences'

subplot(1,2,2)
LocPrecision_hist = 0:0.5:ceil(max(LocPrecision));
histogram(LocPrecision, LocPrecision_hist);
xlabel 'Localization precision (nm)'
ylabel 'Occurences'
xlim([0, ceil(max(LocPrecision))])
title(['Mean: ',num2str(mean(LocPrecision),'%.1f'),' nm - Mode: ',num2str(mode(round5LocPrecision)),' nm - Median: ',num2str(median(round5LocPrecision)),' nm'...
    , ' - STD: ',num2str(std(LocPrecision),'%.1f'),' nm']);

if AnalyseTracks == 1
    Durations = ParticleInfo(:,5);
    figure('name','Emission durations','Units','normalized', 'Position', [0.1 0.1 0.8 0.5],'color','white');
    subplot(1,2,1)
    histogram(Durations);
    xlabel 'Duration (Frames)'
    ylabel 'Occurences'
    title(['Mean: ',num2str(mean(Durations),'%.1f'),' frames - Mode: ',num2str(mode(Durations)),' frames - Median: ',num2str(median(Durations)),' frames'...
        , ' - STD: ',num2str(std(Durations),'%.1f'),' frames']);
    
    subplot(1,2,2)
    plot(ADCs, Durations, '+');
    xlabel 'ADCs'
    ylabel 'Durations (frame)'
    grid on
    
end




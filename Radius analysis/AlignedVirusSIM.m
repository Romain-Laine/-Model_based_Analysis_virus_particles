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

% DATA = 'Simulated';
DATA = 'Measured';

%------------------------------------------------------------------------------------------------------------------------%
ProtStructure = 'Shell';
% ProtStructure = 'Ring';

LinkerStructure = 'Shell';
% LinkerStructure = 'Ring';

radVP = 80;             % nm
ThicknessVP = 30;       % nm

Linker_size = 10;       % nm
Linker_thickness = 2;   % nm
LocError = 15;          % nm

% Total number of localizations
N_loc = 10000000;    % number of localizations (assuming 1 fluorophore / protein)

%------------------------------------------------------------------------------------------------------------------------%
% Fitting parameters
Weighted = 'on';
% Weighted = 'off';

% Number of localizations to generate for the minimization
N_loc_FMS_init = 1000000;
N_iteration_FMS_init = 20; % Maximum number of iterations

N_loc_FMS      = 1000000;   % the number of localization used for the simulation, 10,000,000 seems like a good number
N_iteration_FMS      = 20; % Maximum number of iterations

% Tolerances (these do not seem to do anything !) common to both first and
% second FMS
Tolerance = 0.01;
Chi2_tolerance = 1e-5;

% Initial fitting for better initial guess
N_loc_NLLS_init = 1000000;
N_iteration_NLLS_init = 15;
lambda_init = 0.3;
% Nudging steps for numerical derivative calculation
r_step_init = 1;   % nm
t_step_init = 4;   % nm

N_loc_NLLS      = 10000000;   % the number of localization used for the simulation, 10,000,000 seems like a good number
N_iteration_NLLS      = 10; % Maximum number of iterations
lambda = 0.3;    % damping factor (because it otherwise tends to overshoot)
% Nudging steps for numerical derivative calculation
r_step = 0.5;   % nm
t_step = 2;   % nm

% Minimum number of iteration (common to both first and second NLLS)
Min_N_iteration = 5;
Criterium_dr = 0.01; % Average changes in dr over the last Critera_step number of iterative steps (0.01 = 1% of change)
Criterium_dt = 0.01; % Average changes in dr over the last Critera_step number of iterative steps (0.01 = 1% of change)
Criterium_step = 5;

%------------------------------------------------------------------------------------------------------------------------%
% Pixel size for density display (display only)
PixelSize = 10;
N_display = 10000000; % number of localization to simulate for the display

%------------------------------------------------------------------------------------------------------------------------%
% Standard error and confidence interval
Conf = 0.95; % 95% confidence interval

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
r_hist = linspace(0,Max_r,min(round(N_loc/4),round(Max_r)));
x_hist = linspace(-Max_r,Max_r,min(round(N_loc/4),round(Max_r)));

figure('Color','white','name','Radius histogram');
subplot(1,2,1)
hist(r,r_hist)
xlim([0,Max_r])
xlabel 'Radius (nm)'
ylabel 'Occurences'
title(['Mean radius: ',num2str(mean(r),4),' and STD: ',num2str(std(r),4)])

subplot(1,2,2)
hist(xy(:,1),x_hist)
xlim([-Max_r,Max_r])
xlabel 'x (nm)'
ylabel 'Occurences'

% Display the intensity map
PixelNumber = round(Max_r/PixelSize);
Size = PixelNumber*PixelSize;
SRimage_histogram(xy, PixelNumber,Size);

%% rapidSTORM display
xy_rapidS = xy;
xy_rapidS(:,1) = xy_rapidS(:,1) + 200;
xy_rapidS(:,2) = xy_rapidS(:,2) + 200;

ScaleBarSize = 100; % nm
% SR_imageSize = max(xy_rapidS(:))/PixelSize + 5;
SR_imageSize = 40;
[ SR_image_rS ] = SRdisplay_rapidSTORM( xy_rapidS, PixelSize, SR_imageSize );
ColorImage_rS = Grey2Color(SR_image_rS/max(SR_image_rS(:)));
ColorImage_rS = AddScaleBar( ColorImage_rS, PixelSize, ScaleBarSize);
figure;
imshow(ColorImage_rS)
title(['Scale bar: ',num2str(ScaleBarSize),' nm (',num2str(PixelSize),' nm / pixel)'])
%

disp('----------------------------');


%%

% Linker_size = 30;       % nm
% Linker_thickness = Linker_size/4;   % nm

%------------------------------------------------------------------------------------------------------------------------%
% Calculate the data histogram

% Maximum number of bins permissible
n_bin_max = round(MaxRad); % bins are in unit of nm (no point having smaller bin psacing than 1 nm)
% Number of bins used for the histogramming
% n_bin = min([round(N_loc/4),round(Max_r),n_bin_max]);

% If the data is thresholded at the source then there is no need to use Max_r
n_bin = min([round(N_loc/4),n_bin_max]);

% The bins are chosen to be linearly spaced (is there a wiser choice?)
r_hist = linspace(0,Max_r,n_bin);
N_data = hist(r,r_hist);

% Eliminate the last bin as it may contain all the ouliers in measured
% dataset --> skews results
N_data(end) = [];

disp(['Total number of localizations used for histogramming: ',num2str(sum(N_data))]);
disp(['Number of bins in histogram: ',num2str(n_bin)]);
disp('----------------------------');


%------------------------------------------------------------------------------------------------------------------------%
% Calculate initial values from the histogram itself
% NB: it seems sensible to give both too high or too low values for radius
% and thickness, judging the negative correlation between the 2 parameters

r_init = mean(r);
% t_init = std(r) - LocError/2;
t_init = std(r);

% initialize sim_param variable 
Sim_param = cell(5,1);
Sim_param{1} = ProtStructure;
Sim_param{2} = LinkerStructure;
Sim_param{3} = Linker_size;        % nm
Sim_param{4} = Linker_thickness;   % nm
Sim_param{5} = LocError;           % nm


%% Minimisation using fminsearchs

% ------------------------------------------------------------------------------------------------
tic
disp('First minimization FMS(r,t)...')
disp(['Initial guesses from mean and std: ',num2str([r_init, t_init],4)]);

% Define the anonymous function calculating the reduced Chi2
h_Chi2_fun = @(rt0) CalcChi2( rt0, Sim_param, N_loc_FMS_init, r_hist, N_data);

% Define initial values
rt_init = [r_init t_init];
options = optimset('Display','final','PlotFcns',@optimplotfval,'MaxIter',N_iteration_FMS_init,'TolX',Tolerance,'TolFun',Chi2_tolerance);
SetTol = optimget(options,'TolX');
SetChi2Tol = optimget(options,'TolFun');
% disp(['Defined tolerance: ', num2str(SetTol)]);
% disp(['Defined Chi2 tolerance: ', num2str(SetChi2Tol)]);

% Perform fminsearch
rt_opt = fminsearch(h_Chi2_fun,rt_init,options);

r_init = rt_opt(1);
t_init = rt_opt(2);
disp(['Results after FMS(r,t): ',num2str([r_init, t_init],4)]);
toc
disp('----------------------------');


% Only minimize with thickness

tic
disp('Second minimization FMS(t)...')
disp(['Initial guesses after FMS(r,t): ',num2str([r_init, t_init],4)]);

% Define the anonymous function calculating the reduced Chi2
h_Chi2_fun = @(t0) CalcChi2_singleParam( r_init, t0, Sim_param, N_loc_FMS, r_hist, N_data);

% Initial values are the results from the previous FMS
options = optimset('Display','off','PlotFcns',@optimplotfval,'MaxIter',N_iteration_FMS,'TolX',Tolerance,'TolFun',Chi2_tolerance);
[t_opt,fval] = fminsearch(h_Chi2_fun,t_init,options);

t_init = t_opt; % Change t_init to the value found by FMS(t)
disp(['Fitted radius and thickness after FMS(t):      ',num2str([r_init, t_init],4)]);
toc
disp('----------------------------');


% NLLS fitting
% ------------------------------------------------------------------------------------------------
tic
disp('First NLLS minimization...');
disp(['Initial guess radius and thickness:      ',num2str([r_init, t_init],4)]);

% Initial guesses as found by initial FMS
r_c = r_init;
t_c = t_init;

% Weighting matrix
if strcmp(Weighted, 'on')
    W = N_data;
    W(W == 0) = 1; % Get rid of zero for non-infinite weighting
    W = diag(1./W);
elseif strcmp(Weighted, 'off')
    W = diag(ones(size(N_data,2),1));  % create identity matrix
end

% profile on
h_wait = waitbar(0,'Please wait...');

% Initialize variables to 100% (=1)
av_dr = 1;
av_dt = 1;

% Initialize count to zero
Count = 0;

% The lists will ever be as big as the total maximum number of iterations
r_list = zeros(N_iteration_NLLS_init + N_iteration_NLLS,1);
t_list = zeros(N_iteration_NLLS_init + N_iteration_NLLS,1);
Chi2_list = zeros(N_iteration_NLLS_init + N_iteration_NLLS,1);
dr_list = zeros(N_iteration_NLLS_init + N_iteration_NLLS,1);
dt_list = zeros(N_iteration_NLLS_init + N_iteration_NLLS,1);

% Calculate the degree of freedom for CHi2 calculation
n_freedom = size(N_data,2) - 2 - 1;

% Initial plot without fitted curve
h_histo = figure('Color','white','name','Radius histogram + fitted curve'); % for display during iteration
plot(r_hist(1:end-1),N_data,'+')
xlabel 'Radius (nm)'
ylabel 'Occurences'
grid on

% Make sure it does the number of iteration at least the minimum number set
% by user
if N_iteration_NLLS_init < Min_N_iteration
    N_iteration_NLLS_init = Min_N_iteration;
end

while ( (abs(av_dr) > Criterium_dr || abs(av_dt) > Criterium_dt || Count < Min_N_iteration) && Count < N_iteration_NLLS_init)
    
    % This uses one single wait bar for both NLLS fitting
    waitbar(Count/(N_iteration_NLLS_init + N_iteration_NLLS + 1),h_wait);
    
    % Central position
    N_sim = VirusSim(Sim_param, r_c, t_c, N_loc_NLLS_init, r_hist);
    N_sim_c = N_sim*sum(N_data)/sum(N_sim);
    
    % Calculate the reduced Chi2
    Chi2_c = (1/n_freedom)*sum(((N_sim_c - N_data).^2).*(diag(W))');
    
    % Upper bound with respect to the parameter r (radius)
    N_sim = VirusSim(Sim_param, r_c + r_step_init, t_c, N_loc_NLLS_init, r_hist);
    N_sim_rUB = N_sim*sum(N_data)/sum(N_sim);
    
    % Lower bound with respect to the parameter r (radius)
    N_sim = VirusSim(Sim_param, r_c - r_step_init, t_c, N_loc_NLLS_init, r_hist);
    N_sim_rLB = N_sim*sum(N_data)/sum(N_sim);
    
    % Upper bound with respect to the parameter t (thickness)
    N_sim = VirusSim(Sim_param, r_c, t_c + t_step_init, N_loc_NLLS_init, r_hist);
    N_sim_tUB = N_sim*sum(N_data)/sum(N_sim);
    
    % Lower bound with respect to the parameter t (thickness)
    N_sim = VirusSim(Sim_param, r_c, t_c - t_step_init, N_loc_NLLS_init, r_hist);
    N_sim_tLB = N_sim*sum(N_data)/sum(N_sim);
    
    % Calculate the Jacobian from the numerical derivatives
    J = [(N_sim_rUB' - N_sim_rLB')/(2*r_step_init), (N_sim_tUB' - N_sim_tLB')/(2*t_step_init)];
    
    % Calculate the transpose matrix of the Jacobian matrix
    Jt = transpose(J);
    
    % Calculate the matrix difference
    dN = N_data' - N_sim_c';
    
    % Solve the non-linear least square equation
    dB = (Jt*W*J + lambda_init*diag(diag(Jt*W*J)))\(Jt*W*dN);
    
    dr = dB(1)/r_c;  % dr defined in fraction of r_c
    dt = dB(2)/t_c;  % dt defined in fraction of t_c
    
    % Move the optimal fit by the correct amount
    r_c = r_c + dB(1);
    t_c = t_c + dB(2);
    
    if t_c < 0
        t_c = 0;
    end
    
    % Increment and save
    Count = Count + 1;
    r_list(Count) = r_c;
    t_list(Count) = t_c;
    Chi2_list(Count) = Chi2_c;
    
    disp(['Iteration: ',num2str(Count),' (r,t) = ',num2str([r_c,t_c],4), ' -- Chi2 = ',num2str(Chi2_c)]);
    
    dr_list(Count) = dr;
    dt_list(Count) = dt;
    
    if Count > Criterium_step+1
        av_dr = mean(abs(dr_list(Count - Criterium_step:Count)));
        av_dt = mean(abs(dt_list(Count - Criterium_step:Count)));
    else
        av_dr = mean(abs(dr_list(1:Count)));
        av_dt = mean(abs(dt_list(1:Count)));
    end
    
    % Plot the results from the iteration step
    figure(h_histo);
    plot(r_hist(1:end-1),N_data,'+',r_hist(1:end-1),N_sim_c,'r');
    xlabel 'Radius (nm)'
    ylabel 'Occurences'
    grid on
    title(['Iteration: ',num2str(Count),' (r,t) = ',num2str([r_c,t_c],4)])
    
    %     myFrame = getframe(gcf);
    %     myIm = myFrame.cdata;
    %     imwrite(myIm,[Save_path,'image1_',int2str(Count),'.png']);
    %
end

% profile viewer

% Reassign the actual number of iteration in the initial step of NLLS fitting
N_iteration_NLLS_init = Count;
% Close the figure used for wait bar
disp(['Results after first round of NLLS: ',num2str([r_c, t_c],4)]);
toc
disp('----------------------------');

% Second round of NLLS
% ------------------------------------------------------------------------------------------------
tic
disp('Second NLLS minimization...');
disp(['Initial guesses: ',num2str([r_c, t_c],4)]);

% Initialize the radius and thickness variability variables
av_dr = 1;
av_dt = 1;

figure(h_histo); % for display during iteration
plot(r_hist(1:end-1),N_data,'+')
xlabel 'Radius (nm)'
ylabel 'Occurences'
grid on

% Make sure it does the number of iteration at least the minimum number set
% by user
if N_iteration_NLLS < Min_N_iteration
    N_iteration_NLLS = Min_N_iteration;
end


while ( (abs(av_dr) > Criterium_dr || abs(av_dt) > Criterium_dt || Count < N_iteration_NLLS_init + Min_N_iteration) && Count < N_iteration_NLLS_init + N_iteration_NLLS)
    
    waitbar(Count/(N_iteration_NLLS_init + N_iteration_NLLS + 1),h_wait);
    
    % Central position
    N_sim = VirusSim(Sim_param, r_c, t_c, N_loc_NLLS, r_hist);
    N_sim_c = N_sim*sum(N_data)/sum(N_sim);
    
    % Calculate the reduced Chi2
    Chi2_c = (1/n_freedom)*sum(((N_sim_c - N_data).^2).*(diag(W))');
    
    % Upper bound with respect to the parameter r (radius)
    N_sim = VirusSim(Sim_param, r_c + r_step, t_c, N_loc_NLLS, r_hist);
    N_sim_rUB = N_sim*sum(N_data)/sum(N_sim);
    
    % Lower bound with respect to the parameter r (radius)
    N_sim = VirusSim(Sim_param, r_c - r_step, t_c, N_loc_NLLS, r_hist);
    N_sim_rLB = N_sim*sum(N_data)/sum(N_sim);
    
    % Upper bound with respect to the parameter t (thickness)
    N_sim = VirusSim(Sim_param, r_c, t_c + t_step, N_loc_NLLS, r_hist);
    N_sim_tUB = N_sim*sum(N_data)/sum(N_sim);
    
    % Lower bound with respect to the parameter t (thickness)
    N_sim = VirusSim(Sim_param, r_c, t_c - t_step, N_loc_NLLS, r_hist);
    N_sim_tLB = N_sim*sum(N_data)/sum(N_sim);
    
    % Calculate the Jacobian from the numerical derivatives
    J = [(N_sim_rUB' - N_sim_rLB')/(2*r_step), (N_sim_tUB' - N_sim_tLB')/(2*t_step)];
    
    % Calculate the transpose matrix of the Jacobian matrix
    Jt = transpose(J);
    
    % Calculate the matrix difference
    dN = N_data' - N_sim_c';
    
    % Solve the non-linear least square equation
    dB = (Jt*W*J + lambda*diag(diag(Jt*W*J)))\(Jt*W*dN);
    
    dr = dB(1)/r_c;  % in % for testing
    dt = dB(2)/t_c;  % in % for testing
    
    r_c = r_c + dB(1);
    t_c = t_c + dB(2);
    
    if t_c < 0
        t_c = 0;
    end
    
    
    Count = Count + 1;
    r_list(Count) = r_c;
    t_list(Count) = t_c;
    Chi2_list(Count) = Chi2_c;
    
    disp(['Iteration: ',num2str(Count),' (r,t) = ',num2str([r_c,t_c],4),' -- Chi2 = ',num2str(Chi2_c)]);
    
    dr_list(Count) = dr;
    dt_list(Count) = dt;
    
    if Count > Criterium_step
        av_dr = mean(abs(dr_list(Count - Criterium_step:Count)));
        av_dt = mean(abs(dt_list(Count - Criterium_step:Count)));
    else
        av_dr = mean(abs(dr_list(1:Count)));
        av_dt = mean(abs(dt_list(1:Count)));
    end
    
    % Plot the results from the iteration step
    figure(h_histo);
    plot(r_hist(1:end-1),N_data,'+',r_hist(1:end-1),N_sim_c,'r')
    xlabel 'Radius (nm)'
    ylabel 'Occurences'
    grid on
    title(['Iteration: ',num2str(Count),' (r,t) = ',num2str([r_c,t_c],4)])
    
%     myFrame = getframe(gcf);
%     myIm = myFrame.cdata;
%     imwrite(myIm,[Save_path,'image2_',int2str(Count),'.png']);
    
end

disp(['Count: ',num2str(Count)]);

% profile viewer
close(h_wait);
disp(['Fitted radius and thickness:      ',num2str([r_c, t_c],4)]);
toc


% ---------------------------------------------------------------------------------------------------------
r_opt = r_c;
t_opt = t_c;

r_list(Count+1:end) = [];
t_list(Count+1:end) = [];
Chi2_list(Count+1:end) = [];
dr_list(Count+1:end) = [];
dt_list(Count+1:end) = [];

%

% Display the results of the iteration
figure('Color','white','Position',[100, 100, 1500, 600],'name','Parameters vs. iteration');
subplot(2,3,1)
plot(r_list)
xlabel 'Iteration number'
ylabel 'Radius (nm)'
grid on

subplot(2,3,2)
plot(t_list)
xlabel 'Iteration number'
ylabel('Thickness (nm)')
grid on

subplot(2,3,3)
plot(Chi2_list)
xlabel 'Iteration number'
ylabel('Chi2')
grid on

subplot(2,3,4)
plot(dr_list)
xlabel 'Iteration number'
ylabel('dr')
grid on

subplot(2,3,5)
plot(dt_list)
xlabel 'Iteration number'
ylabel('dt')
grid on

%% Display the results for the optimal parameters
r_opt = 80;
t_opt = 30;

N_opt = VirusSim(Sim_param, r_opt, t_opt, N_display, r_hist);
N_opt = N_opt*sum(N_data)/sum(N_opt);

figure('Color','white','name','Radius histogram - optimal parameters');
plot(r_hist(1:end-1),N_data,'+')
hold on
plot(r_hist(1:end-1),N_opt,'r')
xlabel 'Radius (nm)'
ylabel 'Occurences'
grid on
title(['Fitted radius: ',num2str(r_opt,'%3.2f'), ' nm - Fitted thickness: ',num2str(t_opt,'%3.2f'),' nm'])

%%
%
disp(' ');
disp('--------------------------------------------------');
disp('SUMMARY:');
disp(['Dataset type: ',DATA]);
if strcmp(DATA,'Measured')
    disp(FolderName);
end

disp('Fixed parameters:');
disp(['Linker size: ',num2str(Linker_size),' nm +/- ',num2str(Linker_thickness), ' nm']);
disp(['Localization error: ',num2str(LocError),' nm']);
disp(['Total number of localizations in data: ',num2str(sum(N_data))]);

disp(' ');
disp('RESULTS:')
disp(['Chi^2 = ',num2str(Chi2_c)]);
disp(['Fitted radius: ',num2str(r_opt,'%3.2f'), ' nm - Fitted thickness: ',num2str(t_opt,'%3.2f'),' nm']);


%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Error analysis %%%%%%%%%%%%%%%%%%%%
disp('--------------------------------------------------');
disp(' ');
disp('Calculating standard error on the parameters...');
disp(['Confidence interval: ',num2str(100*Conf),'%']);
N_loc_se = 10000000;

r_step = 0.25; % nm
t_step = 1; % nm

tic
% Upper bound with respect to the parameter r (radius)
N_sim = VirusSim(Sim_param, r_opt + r_step, t_opt, N_loc_se, r_hist);
N_sim_rUB = N_sim*sum(N_data)/sum(N_sim);

% Lower bound with respect to the parameter r (radius)
N_sim = VirusSim(Sim_param, r_opt - r_step, t_opt, N_loc_se, r_hist);
N_sim_rLB = N_sim*sum(N_data)/sum(N_sim);

% Upper bound with respect to the parameter t (thickness)
N_sim = VirusSim(Sim_param, r_opt, t_opt + t_step, N_loc_se, r_hist);
N_sim_tUB = N_sim*sum(N_data)/sum(N_sim);

% Lower bound with respect to the parameter t (thickness)
N_sim = VirusSim(Sim_param, r_opt, t_opt - t_step, N_loc_se, r_hist);
N_sim_tLB = N_sim*sum(N_data)/sum(N_sim);


% Weight array
if strcmp(Weighted, 'on')
    W = N_data;
    W(W == 0) = 1;
    W = 1./W;
elseif strcmp(Weighted, 'off')
    W = ones(1,size(N_data,2));
end

% Calculating the A matrix (see Partick Barrie lecture notes on "Error analysis and curve fitting")
A = zeros(2,2);
A(1,1) = sum(W.*((N_sim_rUB - N_sim_rLB)/(2*r_step)).^2);
A(2,2) = sum(W.*((N_sim_tUB - N_sim_tLB)/(2*t_step)).^2);
A(1,2) = sum(W.*((N_sim_tUB - N_sim_tLB)/(2*t_step)).*((N_sim_rUB - N_sim_rLB)/(2*r_step)));
A(2,1) = A(1,2); % symmetrical

% Residual sum of squares (weighted or not)
RSS = sum(W.*(N_data - N_opt).^2);

n = size(N_data,2); % = n_bin -1 because of gettingn rid of the last bin
p = 2; % number of fitted parameters


% Variance-covariance matrix
V = RSS/(n-p)*inv(A); %#ok<MINV>

% Calculating the confidence interval
alpha = 1-Conf;
ConfInt_r = tinv(1-alpha/2,n-p)*sqrt(V(1,1));
ConfInt_t = tinv(1-alpha/2,n-p)*sqrt(V(2,2));
toc

disp('Standard errors:')
disp(['dr = ',num2str(sqrt(V(1,1)),'%1.3f'), ' nm        dt = ',num2str(sqrt(V(2,2)),'%1.2f'),' nm']);
disp(['Covariance = ',num2str(V(2,1),'%1.3f')]);
disp('Confidence interval:')
disp(['dr = ',num2str(ConfInt_r,'%1.3f'), ' nm        dt = ',num2str(ConfInt_t,'%1.2f'),' nm']);

% Displaying the confidence ellipse (using error_ellipse downloaded from MathWorks)
figure('Color','white','name','Confidence ellipse');
error_ellipse(V,[r_opt,t_opt],'conf',Conf); % see options of error_ellipse for details
hold on 
plot(r_opt,t_opt,'+r')
xlabel 'Radius (nm)'
ylabel 'Thickness (nm)'
grid on
if strcmp(DATA,'Simulated')
    hold on
    plot(radVP,ThicknessVP,'og')
end
axis equal
title([num2str(Conf*100),'% confidence ellipse']);

disp(' ');
disp(['Fitted radius: ',num2str(r_opt,'%3.2f'), ' +/- ',num2str(ConfInt_r,'%3.2f'),' nm - Fitted thickness: ',num2str(t_opt,'%3.2f'),' +/- ',num2str(ConfInt_t,'%3.2f'),' nm']);

% Display how long the whole thing took

dt = toc(tStart);
disp('--------------');
disp(['Elapsed time: ',num2str(dt/60,'%2.1f'),' minute(s)']);

%
beep
pause(0.5);
beep
pause(0.5);
beep
pause(0.5);

% This code allows for the analysis of the assymetry data.
% It does the following:
% - align all the localizations obtained from the mTurquoise localization data.
% - compute the distance of the Centre of Mass of the particle from the aligned capsid (mTurquoise)
% - by fitting the model of a shell (flexibility = thickness), extract the model parameter
% - calculate the 95% confidence interval on the fitted parameters

clear all
close all
clc

Cut_off = 250;   % nm : maximum distance between the 2 centroids
Max_rad = 50;    % nm
Min_ADC = 3400;  % 3,400 ADC represents a cut-off at 10 nm for old settings
LocError = 7.6;

N_sim = 50000;   % seems good based on tests

Conf = 0.95;   % 0.95 = 95% confidence interval

DATA = 'Measured';
% DATA = 'Simulated';

disp(['Data type: ',DATA]);
if strcmp(DATA,'Measured')
    DefaultPath = 'E:\STORM data\2014_07_17 Asymmetry analysis\';
    
    % Load the datasets - Centre of mass from Hough transforms
    [FileName,PathName] = uigetfile('*.txt','Select an cHough centre file',DefaultPath);
    disp(['Loading ',PathName,FileName]);
    XY_Centroids = dlmread([PathName, FileName], '\t',1,0);
    n_Centroids = size(XY_Centroids,1);
    disp(['Number of centroids: ',num2str(n_Centroids)]);
    
    % Load the datasets - Centre of mass from mTurquoise dataset (= capsid)
    [FileName,PathName] = uigetfile('*.txt','Select a capsid centroid file',PathName);
    disp(['Loading ',PathName,FileName]);
    LocInfo = Read_LocFile( [PathName, FileName] );
    
    % Save all localization coordinates into a XY vector
    XY_mTu_loc = LocInfo(:,1:2);
    ADC_mTu = LocInfo(:,4);
    XY_mTu_loc(ADC_mTu < Min_ADC,:) = [];
    
    % Total number of capsid localization
    n_Loc = size(XY_mTu_loc,1);
    
    disp(['Number of capsid localizations: ',num2str(n_Loc)]);
    
    
    %% Display on an SR image-like
    
    figure('Color','white','name','Constellation display');
    plot(XY_mTu_loc(:,1),XY_mTu_loc(:,2),'+r');
    hold on
    plot(XY_Centroids(:,1),XY_Centroids(:,2),'+b');
    xlabel 'x (nm)'
    ylabel 'y (nm)'
    axis equal
    legend('Capsids','cHough centres')
    
    %% Sort out the ones that within reasonnable distance of each other
    
    Im_Centroids_associated = zeros(n_Loc,2);
    D_min = zeros(n_Loc,2);
    
    for i = 1:n_Loc
        D = sqrt((XY_Centroids(:,1) - XY_mTu_loc(i,1)).^2 + (XY_Centroids(:,2) - XY_mTu_loc(i,2)).^2);
        Im_Centroids_associated(i,:) = XY_Centroids(D == min(D),:);
        D_min(i) = min(D);
    end
    
    Coord = Im_Centroids_associated - XY_mTu_loc;
    D = sqrt(Coord(:,1).^2 + Coord(:,2).^2);
    
    % Sorted coordinates to work with
    x = Coord(D < Cut_off,1);
    y = Coord(D < Cut_off,2);
    
    % Get rid of points that are clearly out
    
    x_init = mean(x);
    y_init = mean(y);
    r = sqrt((x-x_init).^2 + (y-y_init).^2);
    
    x(r > Max_rad) = [];
    y(r > Max_rad) = [];
end

% Create simulated data -----------------------------------
if strcmp(DATA,'Simulated')
    
    r_sim = 20;
    N_sim_data = 5000;
    x0 = 40;
    y0 = -11;
    
    disp(['x_0 = ',num2str(x0)]);
    disp(['y_0 = ',num2str(y0)]);
    disp(['r_sim = ',num2str(r_sim)]);
    
    % Generate a dataset with Thickness = 0
    xy_sim = MC_Sim_3DShell(N_sim_data,r_sim,0) + MC_Sim_2DLocError(N_sim_data,LocError);
    x = x0 + xy_sim(:,1);
    y = y0 + xy_sim(:,2);
end


%% Define initial guesses
r = sqrt((x - mean(x)).^2 + (y - mean(y)).^2);

% Other alternative are using mode --> probably more accurate but after testing do not
% give stable results.

% Get rid of the ones that are too far away (bad localization etc.) Farther
% that 3* standard deviation is removed
x(r > 3*std(r)) = [];
y(r > 3*std(r)) = [];

x_init = mean(x);
y_init = mean(y);
r = sqrt((x-x_init).^2 + (y-y_init).^2);
r_init = mean(r);

disp(' ');
disp('------------------------------');
disp(['Number of points for analysis: ',num2str(size(x,1))]);
disp('------------------------------');

h1 = Display_Asym_data(x,y,x_init,y_init,r_init);
set(h1,'Name','Localizations (Initial guesses)');

disp('Initial guesses:');
disp(['x_init = ',num2str(x_init)]);
disp(['y_init = ',num2str(y_init)]);
disp(['r_init = ',num2str(r_init)]);


%% Exhaustive search method - first round ---------------------------------------------------------------------------

disp('------------------------------');
disp('Optimization (round I)...');

% Define extension of the parameter search around the initial guesses
x_ext = 20; % nm
y_ext = 20; % nm
r_ext = 10; % nm

% N_xy defines the number of function estimation (number of parameters to test)
N_xy = 20;  % choosing N = 2x extension means a sampling every 1nm
N_r = 10;

% Define the parameter search space
x0 = linspace(-1,1,N_xy)*x_ext + x_init;
y0 = linspace(-1,1,N_xy)*y_ext + y_init;
r0 = linspace(-1,1,N_r)*r_ext + r_init;
r_max = r_init + r_ext;

disp(['x: from ',num2str(min(x0)),' to ',num2str(max(x0)),' with dx = ',num2str(x0(2)-x0(1))]);
disp(['y: from ',num2str(min(y0)),' to ',num2str(max(y0)),' with dx = ',num2str(y0(2)-y0(1))]);
disp(['r: from ',num2str(max(min(r0),0)),' to ',num2str(max(r0)),' with dx = ',num2str(r0(2)-r0(1))]);

% % The extension of the parameter search for r is between 0 and r_max
% r_max = 30;
% r0 = linspace(0,1,N_r + 1)*r_max;
% r0(1) = []; % do not attempt on r = 0

% Define the histogramming process
n_rbin = ceil(size(x,1)/2);
r_hist = linspace(0,r_max,n_rbin);

n_tbin = 5;  % By using 9 bins the quadrant is divided in 8 equal portions
Theta_hist = linspace(-1,1,n_tbin);

SimParam = [N_sim, LocError];

% Perform the full search
[~, OptParam] = FullSearch_Chi2( x0, y0, r0, x, y, SimParam, r_hist, Theta_hist );

x_init = OptParam(1,1);
y_init = OptParam(1,2);
r_init = OptParam(1,3);

%% Exhaustive search method - Second round ---------------------------------------------------------------------------

disp('------------------------------');
disp('Optimization (round II)...');

% Define extension of the parameter search around the initial guesses
x_ext = 10; % nm
y_ext = 10; % nm
r_ext = 5; % nm

% N_xy defines the number of function estimation (number of parameters to test)
N_xy = 20;  % choosing N = 2x extension means a sampling every 1nm
N_r = 10;

% Define the parameter search space
x0 = linspace(-1,1,N_xy)*x_ext + x_init;
y0 = linspace(-1,1,N_xy)*y_ext + y_init;
r0 = linspace(-1,1,N_r)*r_ext + r_init;
r_max = r_init + r_ext;

disp(['x: from ',num2str(min(x0)),' to ',num2str(max(x0)),' with dx = ',num2str(x0(2)-x0(1))]);
disp(['y: from ',num2str(min(y0)),' to ',num2str(max(y0)),' with dx = ',num2str(y0(2)-y0(1))]);
disp(['r: from ',num2str(max(min(r0),0)),' to ',num2str(max(r0)),' with dx = ',num2str(r0(2)-r0(1))]);

% % The extension of the parameter search for r is between 0 and r_max
% r_max = 30;
% r0 = linspace(0,1,N_r + 1)*r_max;
% r0(1) = []; % do not attempt on r = 0

% Define the histogramming process
n_rbin = ceil(size(x,1)/2);
r_hist = linspace(0,r_max,n_rbin);

n_tbin = 5;  % By using 9 bins the quadrant is divided in 8 equal portions
Theta_hist = linspace(-1,1,n_tbin);

SimParam = [N_sim, LocError];

% Perform the full search
[~, OptParam] = FullSearch_Chi2( x0, y0, r0, x, y, SimParam, r_hist, Theta_hist );

x_init = OptParam(1,1);
y_init = OptParam(1,2);
r_init = OptParam(1,3);


%% Exhaustive search method - Third round ---------------------------------------------------------------------------

disp('------------------------------');
disp('Optimization (round III)...');
% 
% x_init = 13.6334;
% y_init = -36.8198;
% r_init = 23.1117;

% Define extension of the parameter search around the initial guesses
x_ext = 5; % nm
y_ext = 5; % nm
r_ext = 10; % nm

% N_xy defines the number of function estimation (number of parameters to test)
N_xy = 20;  % choosing N = 2x extension means a sampling every 1nm
N_r = 40;

% Define the parameter search space
x0 = linspace(-1,1,N_xy)*x_ext + x_init;
y0 = linspace(-1,1,N_xy)*y_ext + y_init;
r0 = linspace(-1,1,N_r)*r_ext + r_init;
r_max = r_init + r_ext;

disp(['x: from ',num2str(min(x0)),' to ',num2str(max(x0)),' with dx = ',num2str(x0(2)-x0(1))]);
disp(['y: from ',num2str(min(y0)),' to ',num2str(max(y0)),' with dx = ',num2str(y0(2)-y0(1))]);
disp(['r: from ',num2str(max(min(r0),0)),' to ',num2str(max(r0)),' with dx = ',num2str(r0(2)-r0(1))]);
disp(' ');
disp('------------------------------');
disp(['Number of points for analysis: ',num2str(size(x,1))]);
disp('------------------------------');

% % The extension of the parameter search for r is between 0 and r_max
% r_max = 30;
% r0 = linspace(0,1,N_r + 1)*r_max;
% r0(1) = []; % do not attempt on r = 0

% Define the histogramming process
n_rbin = ceil(size(x,1)/2);
r_hist = linspace(0,r_max,n_rbin);

n_tbin = 5;  % By using 9 bins the quadrant is divided in 8 equal portions
Theta_hist = linspace(-1,1,n_tbin);

SimParam = [N_sim, LocError];

% Perform the full search
[K, OptParam] = FullSearch_Chi2( x0, y0, r0, x, y, SimParam, r_hist, Theta_hist );

% --------------------------------------------------------------------------------------------------
% Optimal parameters
x0_opt = OptParam(1,1);
y0_opt = OptParam(1,2);
r0_opt = OptParam(1,3);

% Indices of the optimal parameters
i_opt = OptParam(2,1);
j_opt = OptParam(2,2);
k_opt = OptParam(2,3);

% Normalize the Chi2 for Confidence interval analysis

K2_min = min(K(:));
n_bin = n_rbin + n_tbin - 1 ; % total number of data points
% The last bin of the Theta histogram will never contribute as
% always equal to zero
p = 3;              % number of fitted parameters

Max_diff = finv(Conf,p,2*n_rbin-p)*p*K2_min/(2*n_rbin-p);          % See Dr. P. Barrie lecture notes

K_x0 = K(:,j_opt,k_opt);
K_y0 = K(i_opt,:,k_opt);
K_r0 = squeeze(K(i_opt,j_opt,:));

% Normalize the K^2
norm_K_x0 = (K_x0 - K2_min)/Max_diff;
norm_K_y0 = (K_y0 - K2_min)/Max_diff;
norm_K_r0 = (K_r0 - K2_min)/Max_diff;


figure('Color','white','name','Normalised Chi^2');
subplot(1,3,1)
plot(x0,norm_K_x0,x0,ones(1,size(x0,2)),'r')
xlabel 'x0 (nm)'
ylabel 'Norm. Chi^2'
axis tight
ylim([-0.1 1.2])

subplot(1,3,2)
plot(y0,norm_K_y0,y0,ones(1,size(y0,2)),'r')
xlabel 'y0 (nm)'
ylabel 'Norm. Chi^2'
axis tight
ylim([-0.1 1.2])

subplot(1,3,3)
plot(r0(r0>0),norm_K_r0,r0,ones(1,size(r0,2)),'r')
xlabel 'r0 (nm)'
ylabel 'Norm. Chi^2'
axis tight
ylim([-0.1 1.2])


norm_Kxy = (K(:,:,k_opt)-K2_min)/Max_diff;
figure('Color','white','name','Chi^2 space image at r = r_opt');
imshow(norm_Kxy,[]);

% Finding the bounds CI @ 95%
[ x0_min, x0_max ] = Find_bounds(x0, norm_K_x0 );
[ y0_min, y0_max ] = Find_bounds(y0, norm_K_y0'); % The " ' " is needed to give the correct input
[ r0_min, r0_max ] = Find_bounds(r0(r0>0), norm_K_r0 );

disp(' ');
disp('------------------------------');
disp(['Confidence bounds (',num2str(100*Conf),'%):']);
disp(['x0: ', num2str([ x0_min, x0_max ])])
disp(['y0: ', num2str([ y0_min, y0_max ])])
disp(['r0: ', num2str([ r0_min, r0_max ])])
disp('------------------------------');


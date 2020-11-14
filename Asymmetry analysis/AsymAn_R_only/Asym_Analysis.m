% This code allows for the analysis of the assymetry data.
% It does the following:
% - align all the localizations obtained from the mTurquoise localization data.
% - compute the distance of the Centre of Mass of the particle from the aligned capsid (mTurquoise)
% - by fitting the model of a shell (flexibility = thickness), extract the model parameter
% - calculate the 95% confidence interval on the fitted parameters

clear all
close all
clc

% This needs to be here !
D_cut_off = 250;   % nm : maximum distance between the 2 centroids
r_cut_off = 100;

N_sim = 1000000;   % seems good based on tests

Conf = 0.95;   % 0.95 = 95% confidence interval

DATA = 'Measured';
% DATA = 'Simulated';

prompt = {'Enter Min ADC:','Enter number of datasets to load:'};
dlg_title = 'Input';
num_lines = 1;
def = {'3400','4'};
answer = inputdlg(prompt,dlg_title,num_lines,def);
Min_ADC = str2double(answer{1});
N_dataset = str2double(answer{2});

x_lim = zeros(1,4);

disp(['Data type: ',DATA]);
if strcmp(DATA,'Measured')
    DefaultPath = 'E:\STORM data\2014_07_17 Asymmetry analysis\';
    x = zeros(0,0);
    y = zeros(0,0);
    
    h1 = figure('Color','white','name','Constellation display','units','normalized','position',[0.01 0.1 0.98 0.7]);
    
    for j = 1:N_dataset
        % Load the datasets - Centre of mass from Hough transforms
        disp('---------------------');
        [FileName,PathName] = uigetfile('*.txt','Select an cHough centre file',DefaultPath);
        DefaultPath = PathName;
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
        
        subplot(2,N_dataset,j);
        plot(XY_mTu_loc(:,1),XY_mTu_loc(:,2),'+r');
        hold on
        plot(XY_Centroids(:,1),XY_Centroids(:,2),'+b');
        xlabel 'x (nm)'
        ylabel 'y (nm)'
        axis equal
        legend('Capsids','cHough centres')
        xlim([0 41000]);
        ylim([0 41000]);
        
        %% Sort out the ones that within reasonnable distance of each other
        
        Im_Centroids_associated = zeros(n_Loc,2);
        
        % Find the nearest corresponding virus
        for i = 1:n_Loc
            D = sqrt((XY_Centroids(:,1) - XY_mTu_loc(i,1)).^2 + (XY_Centroids(:,2) - XY_mTu_loc(i,2)).^2);
            Im_Centroids_associated(i,:) = XY_Centroids(D == min(D),:);
        end
        
        % Centre on capsid centroid
        xi = Im_Centroids_associated(:,1) - XY_mTu_loc(:,1);
        yi = Im_Centroids_associated(:,2) - XY_mTu_loc(:,2);
        D = sqrt(xi.^2 + yi.^2);
        
        xi(D > D_cut_off) = [];
        yi(D > D_cut_off) = [];
        
        % Centre on the centre of mass
        xi = xi - mean(xi);
        yi = yi - mean(yi);
        
        % Remove the ones that are too far from the Centre of mass
        r = sqrt(xi.^2 + yi.^2);
        xi(r > r_cut_off) = [];
        yi(r > r_cut_off) = [];
        
        % Re-centre the dataset after removing the outliers
        xi = xi - mean(xi);
        yi = yi - mean(yi);
        
        % Remove the ones that are too far from the Centre of mass
        r = sqrt(xi.^2 + yi.^2);
        SD_cutoff = mean(r) + 3*std(r);
        % Remove the datapoints that lie further than 3* STD
        xi(r > SD_cutoff) = [];
        yi(r > SD_cutoff) = [];
        
        disp(['3xSD cut-off: ',num2str(SD_cutoff),' nm']);
        
        %         Re-centre again after having removed the outliers
        xi = xi - mean(xi);
        yi = yi - mean(yi);
        
        disp(['Number of particles useable: ',num2str(size(xi,1))]);
        
        % For display
        x_lim(i) = 1.1*max(cat(1,abs(xi), abs(yi)));
        
        subplot(2,N_dataset,N_dataset+j);
        plot(xi,yi,'+');
        hold on
        plot(0,0,'o')
        axis equal
        grid on
        xlabel 'x (nm)'
        ylabel 'y (nm)'
        title([num2str(size(xi,1)),' particles'])
        
        % Sorted coordinates to work with,  centre on the centre of mass
        x = cat(1, x, xi);
        y = cat(1, y, yi);
        
    end
    
    for j = 1:N_dataset
        subplot(2,N_dataset,N_dataset+j);
        xlim([-max(x_lim) max(x_lim)]);
        ylim([-max(x_lim) max(x_lim)]);
    end
    
    
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


%% -------------------

% Initialize dataset
r = sqrt(x.^2 + y.^2);

SD_cutoff = mean(r) + 2*std(r);
% Remove the datapoints that lie further than 2* STD
x(r > SD_cutoff) = [];
y(r > SD_cutoff) = [];

% Re-centre again after having removed the outliers
x = x - mean(x);
y = y - mean(y);
% Initialize dataset
r = sqrt(x.^2 + y.^2);

h = Display_Asym_data(x,y,0,0,mean(r));
set(h,'name','Localizations');

% Define histogram
r_hist = 0:2:floor(SD_cutoff);
% Every 2 nm between 0 and 3 SD_cutoff

n_r = histc(r,r_hist);

disp(' ');
disp('------------------------------');
disp(['3xSD cut-off: ',num2str(SD_cutoff),' nm']);
disp(['Total number of points for analysis: ',num2str(size(x,1))]);
disp('------------------------------');
% % Define extension of the parameter search around the initial guesses
% r_ext = 10; % nm
% % N_r defines the number of function estimation (number of parameters to test)
% N_r = 200;
% % Define the parameter search space
% r_span = linspace(-1,1,N_r)*r_ext + r_init;

r_span = 0:0.1:30;
disp(['r: from ',num2str(max(min(r_span),0)),' to ',num2str(max(r_span)),' with dr = ',num2str(r_span(2)-r_span(1))]);

SimParam{1} = N_sim;
SimParam{2} = r_span;

%% Exhaustive search method ---------------------------------------------------------------------------

% LocError = [6 6.5 7 7.5 8 8.5 9 9.5 10];
LocError = 8;

Results = zeros(0,5);

for i = 1:size(LocError,2)
    SimParam{3} = LocError(i);
    disp(' ');
    disp('------------------------------');
    disp('Optimization ...');
    disp(['Loc. error = ',num2str(LocError(i)),' nm']);
    
    % Perform the full search
    K_r = FullSearch_Chi2(r_hist, n_r, SimParam);
    
    % --------------------------------------------------------------------------------------------------
    % Optimal parameter and indice of the optimal parameter
    [minK, k_opt] = min(K_r);
    r_opt = r_span(k_opt);
    
    % Normalize the Chi2 for Confidence interval analysis
    
    p = 1;              % number of fitted parameters
    n_rbin = size(r_hist,2);
    Max_diff = finv(Conf,p,2*n_rbin-p)*p*minK/(2*n_rbin-p);          % See Dr. P. Barrie lecture notes
    
    % Normalize the K^2
    norm_K_r = (K_r - minK)/Max_diff;
    
    h2 = figure('Color','white','name','Normalised Chi^2');
    plot(r_span(r_span > 0), norm_K_r, r_span, ones(1,size(r_span,2)),'r')
    xlabel 'r0 (nm)'
    ylabel 'Norm. Chi^2'
    axis tight
    ylim([-0.1 1.2])
    
    % Finding the bounds CI @ 95%
    [ r0_min, r0_max ] = Find_bounds(r_span(r_span > 0), norm_K_r );
    disp(['Confidence bounds (',num2str(100*Conf),'%):']);
    disp(['r0: ', num2str([ r0_min, r0_max ])])
    disp('------------------------------');
    
    Results = cat(1,Results,[LocError(i) r_opt r0_min r0_max minK]);
end

disp(' ');
disp('-------- Results ---------');
disp(Results);

beep;
pause(0.5);
beep;
pause(0.5);
beep;


function [ K, OptParam ] = FullSearch_Chi2( x0, y0, r0, x, y, SimParam, r_hist, Theta_hist )
%FULLSEARCH_CHI2 Summary of this function goes here
%   Detailed explanation goes here

N_sim = SimParam(1);
LocError = SimParam(2);

n_sample = size(x,1);
N_xy = size(x0,2);

r0(r0 <= 0) = []; % Do not attempt fitting where r is negative
N_r = size(r0,2); % Important: by taking N_r = size... the reduction of the size of r0 is taken into account
% Therefore N_r might be different from the initially set N_r

% Minimization criterium (Chi^2)
K_r = zeros(N_xy,N_xy,N_r);
K_t = zeros(N_xy,N_xy);

p = 3; % Number of fitted parameters
n_freedom = size(r_hist,2) + size(r_hist,2) - p - 2;
% The last bin of the Theta histogram will never contribute as
% always equal to zero

% For the Theta histogramming
n_tShell = ones(size(Theta_hist,2),1)*n_sample/(size(Theta_hist,2)-1);
n_tShell(end) = 0;
% Weighting can be chosen to be model weighting
w_t = (size(Theta_hist,2)-1)/n_sample;

h_wait = waitbar(0,'Please wait...');
Count = 0;
tic
for i = 1:N_xy
    Count = Count + 1;
    waitbar(Count/(N_xy),h_wait);
    
    for j = 1:N_xy
        % Calculate the radii with the currect center (x0,y0)
        r = sqrt((x-x0(i)).^2 + (y-y0(j)).^2);
        Theta = atan2(y-y0(j),x-x0(i))/pi;
        Theta(Theta == 1) = -1;
        % All angles meaured as Pi are set to -Pi to match the way histc
        % histograms the data
        
        % Calculate the r histogram
        n_r = histc(r,r_hist);
        w_r = 1./n_r;
        w_r(w_r == Inf) = 1;
        
        % Calculate the Theta histogram
        n_t = histc(Theta,Theta_hist);
        
        K_t(i,j) = sum(w_t.*(n_t - n_tShell).^2);
        % The last bin of the Theta histogram will never contribute as
        % always equal to zero
        
        for k = 1:N_r
            % Generate a dataset with Thickness = 0
            xy_sim = MC_Sim_3DShell(N_sim,r0(k),0) + MC_Sim_2DLocError(N_sim,LocError);
            r_sim = sqrt(xy_sim(:,1).^2 + xy_sim(:,2).^2);
            
            n_rShell = histc(r_sim,r_hist);
            n_rShell = n_sample*n_rShell/sum(n_rShell);
            
            K_r(i,j,k) = sum(w_r.*(n_r - n_rShell).^2);

        end
    end
end
close(h_wait);

% Calculate the total Chi2
K = K_r;
for i = k:N_r
    K(:,:,k) = K(:,:,k) + K_t;
end
K = (1/n_freedom)*K;


% Check the results - find the minimim of Chi^2
[~,min_index] = min(K(:));
[i_opt,j_opt,k_opt] = ind2sub(size(K),min_index);

% Save those variable to base workspace
assignin('base', 'K_r', K_r);
assignin('base', 'K_t', K_t);
assignin('base', 'k_opt', k_opt);

disp('RESULTS:');
disp(['Chi^2 = ',num2str(min(K(:)))]);
disp(['x0 = ',num2str(x0(i_opt)),' nm']);
disp(['y0 = ',num2str(y0(j_opt)),' nm']);
disp(['r0 = ',num2str(r0(k_opt)),' nm']);

OptParam = [x0(i_opt) y0(j_opt) r0(k_opt); i_opt j_opt k_opt];
toc

% Display results
h = Display_Asym_data(x,y,x0(i_opt),y0(j_opt),r0(k_opt));
set(h,'Name','Localizations');

N_display = N_sim*10;
r = sqrt((x-x0(i_opt)).^2 + (y-y0(j_opt)).^2);
Theta = atan2(y-y0(j_opt),x-x0(i_opt))/pi;

n_r = histc(r,r_hist);
n_t = histc(Theta,Theta_hist);

xy_sim = MC_Sim_3DShell(N_display,r0(k_opt),0) + MC_Sim_2DLocError(N_display,LocError);
r_sim = sqrt(xy_sim(:,1).^2 + xy_sim(:,2).^2);
n_rShell = histc(r_sim,r_hist);
n_rShell = n_sample*n_rShell/sum(n_rShell);

n_tShell = ones(size(Theta_hist,2),1)*n_sample/(size(Theta_hist,2)-1);
n_tShell(end) = 0;


% Radius histogram
figure('Color','white','name','Results - Radius histograms');
plot(r_hist,n_r,'+b-')
hold on
plot(r_hist,n_rShell,'r-')
xlabel 'radius (nm)'
grid on
legend('Data histogram','Shell distribution histogram','location','NorthWest')
title(['r_0 = ',num2str(r0(k_opt)), ' nm'])

% Theta histogram
figure('Color','white','name','Results - Theta histograms');
plot(Theta_hist,n_t,'+b-')
hold on
plot(Theta_hist,n_tShell,'r-')
xlabel 'Theta (unit of Pi)'
grid on
legend('Data histogram','Shell distribution histogram','location','NorthWest')
title(['r_0 = ',num2str(r0(k_opt)), ' nm'])

% Plot the Chi2
K_x0 = K(:,j_opt,k_opt);
K_y0 = K(i_opt,:,k_opt);
K_r0 = squeeze(K(i_opt,j_opt,:));

figure('Color','white','name','Results - Chi^2');
subplot(1,3,1)
plot(x0,K_x0)
xlabel 'x0 (nm)'
ylabel 'Chi^2'
axis tight

subplot(1,3,2)
plot(y0,K_y0)
xlabel 'y0 (nm)'
ylabel 'Chi^2'
axis tight

subplot(1,3,3)
plot(r0,K_r0)
xlabel 'r0 (nm)'
ylabel 'Chi^2'
axis tight



end


function [ K_r ] = FullSearch_Chi2(r_hist, n_r, SimParam)
%FULLSEARCH_CHI2 Summary of this function goes here
%   Detailed explanation goes here

N_sim = SimParam{1};
r_span = SimParam{2};
LocError = SimParam{3};

n_sample = sum(n_r);

r_span(r_span <= 0) = []; % Do not attempt fitting where r is negative
N_r = size(r_span,2); % Important: by taking N_r = size... the reduction of the size of r_span is taken into account
% Therefore N_r might be different from the initially set N_r

% Minimization criterium (Chi^2)
K_r = zeros(1,N_r);

p = 1; % Number of fitted parameters
n_freedom = size(r_hist,2) - p - 2;
% The last bin of the Theta histogram will never contribute as
% always equal to zero

% Calculate the r histogram
w_r = 1./n_r;
w_r(w_r == Inf) = 1;

h_wait = waitbar(0,'Please wait...');
tic

for k = 1:N_r
    waitbar(k/N_r);
    % Generate a dataset with Thickness = 0
    xy_sim = MC_Sim_3DShell(N_sim,r_span(k),0) + MC_Sim_2DLocError(N_sim,LocError);
    r_sim = sqrt(xy_sim(:,1).^2 + xy_sim(:,2).^2);
    
    n_rShell = histc(r_sim,r_hist);
    n_rShell = n_sample*n_rShell/sum(n_rShell);
    K_r(k) = (1/n_freedom)*sum(w_r.*(n_r - n_rShell).^2);
    
end

close(h_wait);

[minK, k_opt] = min(K_r);

disp('RESULTS:');
disp(['Chi^2 = ',num2str(minK)]);
disp(['r0 = ',num2str(r_span(k_opt)),' nm']);
toc


N_display = N_sim*10;
xy_sim = MC_Sim_3DShell(N_display,r_span(k_opt),0) + MC_Sim_2DLocError(N_display,LocError);
r_sim = sqrt(xy_sim(:,1).^2 + xy_sim(:,2).^2);
n_rShell = histc(r_sim,r_hist);
n_rShell = n_sample*n_rShell/sum(n_rShell);


% Radius histogram
figure('Color','white','name','Results - Radius histograms/Chi^2','units','normalized','position',[0.01 0.1 0.5 0.4]);
subplot(1,2,1)
plot(r_hist,n_r,'+b-')
hold on
plot(r_hist,n_rShell,'r-')
xlabel 'radius (nm)'
grid on
legend('Data histogram','Shell distribution histogram','location','NorthWest')
title(['r_0 = ',num2str(r_span(k_opt)), ' nm'])


% Plot the Chi2
% figure('Color','white','name','Results - Chi^2');
subplot(1,2,2)
plot(r_span,K_r)
xlabel 'r0 (nm)'
ylabel 'Chi^2'
axis tight



end


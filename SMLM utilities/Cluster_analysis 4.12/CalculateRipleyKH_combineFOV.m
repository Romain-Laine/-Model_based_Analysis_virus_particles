function [ K, H, H2, D ] = CalculateRipleyKH_combineFOV( X_all, Y_all, Area_all, Analysis_window, r_step, Fov_all )
%UNTITLED7 Summary of this function goes here
%   Detailed explanation goes here

% Area should already be in nm^2
% Fov is in nm


tic

N_fov = size(X_all,1);
r_hist = (0:r_step:Analysis_window)';
N_hist_combined = zeros(length(r_hist),1);
D_all = zeros(1,N_fov);

N_loc_1 = 0;
N_loc_2 = 0;

for i = 1:N_fov
    X1 = X_all{i,1};
    X2 = X_all{i,2};
    Y1 = Y_all{i,1};
    Y2 = Y_all{i,2};
    
    % Area should already be in nm^2
    [ r_hist, N_hist, ~, n1, ~ ] = CalcRipleyK( X1, Y1, X2, Y2, Fov_all(i,:), Area_all(i), Analysis_window, r_step );
    N_hist_combined = N_hist_combined + N_hist; % here averaging is done afterwards
    
    % Add the localizations from that FOV to the total FOV loc count
    N_loc_2 = N_loc_2 + length(X2);
    
    D_all(i) = length(X2)/Area_all(i);
    N_loc_1 = N_loc_1 + n1;
    
end

disp('------------------------------------------------');
D = mean(D_all);
disp(['Overal density: ',num2str(10^6*D), ' centroids / um^2']);

% N_hist_combined(end) = []; % Get rid of the last one
K = (1/D)*cumsum(N_hist_combined/N_loc_1); % averaging performed here
H = sqrt(K/pi) - r_hist;  % H(r) = L(r) - r

% H2 outputs the avergae number of neighbours
H2 = D*K - D*pi*r_hist.^2;

toc

end


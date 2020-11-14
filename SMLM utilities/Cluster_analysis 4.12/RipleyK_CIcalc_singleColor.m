function [ H_mean, H_std, Envelope, CI_upperbound, CI_lowerbound ] = RipleyK_CIcalc_singleColor( Fov, Density, Analysis_window, r_step, N_repeat, p )
%RIPLEYK_CICALC Summary of this function goes here
%   Detailed explanation goes here

Display = 0;

tic
disp(['Number of simulations used for calculating CI: ', num2str(N_repeat)]);
disp(['Error type 1: ', num2str(100*p), '%']);
disp(['Density used to calculate CI: ',num2str(10^6*Density), ' localizations / um^2']);

av_N_loc = ceil(Density*Fov(1)*Fov(2));
disp(['Number of localizations used to calculate CI: ',num2str(av_N_loc)]);

r_hist = (0:r_step:Analysis_window)';
K_all = zeros(length(r_hist), N_repeat);
H_all = zeros(length(r_hist), N_repeat);

% For Complete Spatial Randomness (CRS) the number of loc in a specifi FOV
% is given by a Poisson distribution of parameter given by the average
% number of loc / fov

N_loc1 = poissrnd(av_N_loc, N_repeat, 1);

for i = 1:N_repeat
    [ X,Y ] = Generate_random_XY( Fov, N_loc1(i) );
    [ ~, ~, K_all(:,i), ~ ] = CalcRipleyK( X, Y, X, Y, Fov, Fov(1)*Fov(2), Analysis_window, r_step );
    H_all(:,i) = sqrt(K_all(:,i)/pi) - r_hist;  % needs to start at r_step H function is H(r) = L(r) - r
end

H_mean = mean(H_all,2);
H_std = std(H_all,0,2);
Envelope = max(H_all,[],2);
[ CI_upperbound, CI_lowerbound ] = CalculateCIbound( H_all, p );

if Display == 1
    figure('Color','white','name','CI calculations: mean L(r)-r ', 'Units','normalized');
    plot(r_hist(1:end-1), H_all)
    xlim([0 Analysis_window])
    xlabel('r (nm)')
    ylabel('Ripley''s L(r)-r function')
    grid on
    
    figure('Color','white','name','CI calculations: mean L(r)-r ', 'Units','normalized');
    plot(r_hist(1:end-1), H_mean)
    hold on
    plot(r_hist(1:end-1), H_std, '--r')
    hold on
    plot(r_hist(1:end-1), -H_std, '--r')
    xlim([0 Analysis_window])
    xlabel('r (nm)')
    ylabel('Ripley''s L(r)-r function')
    grid on
end

toc

end



function [ CI_upperbound, CI_lowerbound ] = CalculateCIbound( H_all, p )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

N_step = size(H_all,1);
N_sim = size(H_all,2);

CI_upperbound = zeros(N_step,1);
CI_lowerbound = zeros(N_step,1);

for i = 1:N_step
    Sorted_H = sort(H_all(i,:));
    CI_upperbound(i) = Sorted_H(round(N_sim*(1-p/2)));
    CI_lowerbound(i) = Sorted_H(round(N_sim*(p/2)));
end



end


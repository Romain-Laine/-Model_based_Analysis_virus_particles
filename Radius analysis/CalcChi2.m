function [ Chi2 ] = CalcChi2( rt, Sim_param, N_loc_sim, r_hist, n_data)
% CalcChi2 calculates the reduced Chi2 for the minimisation.
% The model is the 3D shell model

r = rt(1);
t = rt(2);

% Generate the model using N localizations
n_sim = VirusSim(Sim_param, r, t, N_loc_sim, r_hist);
n_sim = sum(n_data)*n_sim/(sum(n_sim));


% Calculate the Weighting factor for reduced Chi2
W = n_data;
W(W == 0) = 1;

% Calculate the degree of freedom
n_freedom = size(n_data,2) - 2 - 1;

% Calculate Chi2 with the correct degree of freedom 
Chi2 = (1/n_freedom)*sum((n_data - n_sim).^2./W);

end


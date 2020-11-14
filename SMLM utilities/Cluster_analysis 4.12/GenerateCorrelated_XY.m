function [ X_corr, Y_corr ] = GenerateCorrelated_XY( X,Y, Loc_distance, Loc_distance_flexibility, f)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here

N_loc_corr = round(size(X,1)*f/100);

Theta = 2*pi*rand(N_loc_corr,1);
distance_dc = Loc_distance + Loc_distance_flexibility*(rand(N_loc_corr,1)-0.5);

X_corr(1:N_loc_corr) = X(1:N_loc_corr) + distance_dc.*cos(Theta);
Y_corr(1:N_loc_corr) = Y(1:N_loc_corr) + distance_dc.*sin(Theta);
X_corr = X_corr';
Y_corr = Y_corr';

end


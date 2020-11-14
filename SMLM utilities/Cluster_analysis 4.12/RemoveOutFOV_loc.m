function [ Xr, Yr ] = RemoveOutFOV_loc( X,Y, Fov )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

disp('Removing localizations outside of the FOV...');

Xr = X;
Yr = Y;
OutFOV_loc = (X < 0| Y < 0 | X > Fov(1) | Y > Fov(2));   
Xr(OutFOV_loc) = [];
Yr(OutFOV_loc) = [];
disp(['Number of localizations removed: ',num2str(length(X)-length(Xr))]);


end


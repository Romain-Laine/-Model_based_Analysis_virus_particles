function [ X,Y ] = Generate_random_XY( Fov, N_loc )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

Display_info = 0;

X = Fov(1)*rand(N_loc,1);
Y = Fov(2)*rand(N_loc,1);

if Display_info == 1
    disp(['Simulated density: ',num2str(N_loc*1e6/(Fov(1)*Fov(2))),' particle / um^2']);
end

end


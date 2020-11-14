function [ xy ] = MC_Sim_2DCircle(N,R)
% MC_Sim_2DCircle Monte-Carlo simulation of a 2D circle
% N is the number of point to pick on a ring
% R is the radius of the circle


Theta = 2*pi*rand(N,1);
xy = [R*cos(Theta), R*sin(Theta)];


end


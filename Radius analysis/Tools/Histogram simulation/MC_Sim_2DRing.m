function [ xy ] = MC_Sim_2DRing(N, R, T)
% MC_Sim_2DRing Monte-Carlo simulation of a 2D ring
% N is the number of point to pick on a ring
% R is the radius of the ring
% T is the thickness of the ring

Theta = 2*pi*rand(N,1);

if T ~= 0
    R_max = R + T/2;
    R_min = R - T/2;
    r = sqrt((R_max^2 - R_min^2)*rand(N,1) + R_min^2);
else
    r = R;
end

xy = [r.*cos(Theta), r.*sin(Theta)];


end


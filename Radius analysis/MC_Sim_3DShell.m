function [ xy ] = MC_Sim_3DShell(N, R, T)
% MC_Sim_3DShell Monte-Carlo simulation of a 3D shell projected on a 2D
% plane
% N is the number of point to pick on a shell
% R is the radius of the shell
% T is the thickness of the shell

if T > 0
    R_max = R + T/2;
    R_min = R - T/2;
    r = ((R_max^3 - R_min^3)*rand(N,1) + R_min^3).^(1/3);       % radius
else
    r = R;
end

Phi = 2*pi*rand(N,1);                                       % polar
Theta = acos(2*rand(N,1) - 1);                              % azimuth
xy = [r.*sin(Theta).*cos(Phi), r.*sin(Theta).*sin(Phi)];


end


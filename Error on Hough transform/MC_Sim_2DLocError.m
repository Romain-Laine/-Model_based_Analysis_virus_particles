function [ xy ] = MC_Sim_2DLocError(N,LocError)
% MC_Sim_2DLocError Monte-Carlo simulation of a 2D gaussian localization
% error
% N is the number of point to pick on a ring
% LocError is the standard deviation of the Gaussian localization precision
% error

if LocError ~= 0
    Theta = 2*pi*rand(N,1);
    r = sqrt(2)*LocError*randn(N,1); % Distance of localisation error
    % sqrt(2) factor is for adopting the same convention as Thompson 2002
    % formulae and rainSTORM
    xy = [r.*cos(Theta), r.*sin(Theta)];
else
    xy = zeros(N,2);
end


end


function [ X_d, Y_d ] = SimTrafo( X,Y, Trafo_param)
%SIMTRAFO Summary of this function goes here
%   Detailed explanation goes here

Trafo_type = Trafo_param{1};
Defo_centre = Trafo_param{2};
Defo_Amp_lin = Trafo_param{3};
Defo_Amp_quad = Trafo_param{4};
xy_shift = Trafo_param{5};

disp(['Simulating trafo: ',Trafo_type]);

X_d = X - Defo_centre(1);
Y_d = Y - Defo_centre(2);

% apply the centro-symmetric circular transform
[Theta_d,r_d] = cart2pol(X_d,Y_d);
r_d = r_d + Defo_Amp_lin*r_d + Defo_Amp_quad*r_d.^2 ;
[X_d,Y_d] = pol2cart(Theta_d,r_d);

% re-centre the coordinates and apply the translational shift
X_d = X_d + Defo_centre(1) + xy_shift(1);
Y_d = Y_d + Defo_centre(2) + xy_shift(2);



end


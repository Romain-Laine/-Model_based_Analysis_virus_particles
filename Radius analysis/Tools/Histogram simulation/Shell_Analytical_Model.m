function rHisto = Shell_Analytical_Model(r_hist,radVP,ThicknessVP,LocError,Linker_size, Linker_width)
% Shell_Analytical_Model computes the radius distribution of a projected shell.
% This model does not take into account the effect of the linker or the
% localization error. The calculation is based on the cylindrical
% coordinates infinitesimal volume: dV = r*dr*dPhi*dz


tic
LocError2D = LocError; % Using the image method for the localisation "PSF"
                       % we want PSF = exp-(x^2 + y^2)/(2s^2)
                       % Where s = Thompson(2D) localisation error
Sim_Step = 0.5;   % nm


% ----------------------------------------------------------------
x_prot = -(radVP+ThicknessVP/2):Sim_Step:(radVP+ThicknessVP/2);
[X_prot,Y_prot] = meshgrid(x_prot);
Z_prot = real(sqrt((radVP+ThicknessVP/2)^2 - X_prot.^2 - Y_prot.^2)) - real(sqrt((radVP-ThicknessVP/2)^2 - X_prot.^2 - Y_prot.^2));


% ---------------
x_linker = -(Linker_size + Linker_width/2):Sim_Step:(Linker_size + Linker_width/2);
[X_linker,Y_linker] = meshgrid(x_linker);
Z_linker = real(sqrt((Linker_size+Linker_width/2)^2 - X_linker.^2 - Y_linker.^2)) - real(sqrt((Linker_size-Linker_width/2/2)^2 - X_linker.^2 - Y_linker.^2));


% ---------------
x_LocError = -(4*LocError):Sim_Step:(4*LocError);
[X_LocError,Y_LocError] = meshgrid(x_LocError);
Z_LocError = exp(-(X_LocError.^2 + Y_LocError.^2)/(2*LocError2D^2));
           % Probably fine for small Sim_Step
           % Maybe use erf() for large Sim_Step


% Convolution
Z_tot = conv2(Z_LocError,conv2(Z_prot,Z_linker));

%%
n = size(Z_tot,1);
MyRow = Z_tot((n-1)/2+1,(n-1)/2+1:end);

r = 0:Sim_Step:Sim_Step*(n-1)/2;
Shell_model = abs(r).*MyRow;
rHisto = interp1(r,Shell_model,r_hist);
rHisto = rHisto/sum(rHisto);

disp('Analytical model:');
toc

% assignin('base','modelr',r);
% assignin('base','modelShell_model',Shell_model);

end


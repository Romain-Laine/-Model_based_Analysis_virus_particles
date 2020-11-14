clc
clear all
close all


Step = 1; % nm

radVP = 80;
Thickness = 30;

Linker_size = 10;
Linker_flex = 2;

Sigma = 15;

Max_xy = radVP + Thickness/2 + Linker_size + Linker_flex/2 + Sigma;
gv = -Max_xy:Step:Max_xy;  % nm
[X,Y,Z] = meshgrid(gv);

R = sqrt(X.^2 + Y.^2 + Z.^2);
Sphere = ones(size(R));
Sphere(R > radVP + Thickness/2) = 0;
Sphere(R < radVP - Thickness/2) = 0;

Sphere = smooth3(Sphere,'gaussian',11,4);

Sphere_cut = Sphere;
% Sphere_cut((X >= 0) & (Y <= 0) & (Z >= 0)) = 0;   % for a quarter of the
% shell open
Sphere_cut((Y <= 0) & (Z >= 0)) = 0;

figure('Color','white');

disp('Isosurface...');
tic
isoval = 0.5;
hisoxy = patch(isosurface(Sphere_cut,isoval),...
	'FaceColor','blue',...
	'EdgeColor','none');
isonormals(Sphere_cut,hisoxy)
toc

set(gcf,'Renderer','zbuffer'); lighting phong
set(hisoxy,'SpecularColorReflectance',0.2,'SpecularExponent',50)

lightangle(30,30);
view(24, 12)
axis off
grid off
axis equal


%%

gv_small = -(Linker_size + Linker_flex/2):Step:(Linker_size + Linker_flex/2);
[X_small,Y_small,Z_small] = meshgrid(gv_small);
R_small = sqrt(X_small.^2 + Y_small.^2 + Z_small.^2);

SmallSphere = ones(size(R_small));
SmallSphere (R_small > Linker_size + Linker_flex/2) = 0;
SmallSphere (R_small < Linker_size - Linker_flex/2) = 0;

disp('Convolution...');
tic
Sphere_conv = convn(Sphere,SmallSphere,'same');
toc

Sphere_conv_cut = Sphere_conv;
Sphere_conv_cut((X >= 0) & (Y <= 0) & (Z >= 0)) = 0;
Sphere_conv_cut = smooth3(Sphere_conv_cut);
figure('Color','white');
disp('Isosurface...');
tic
isoval = 0.5;
hisoxy = patch(isosurface(Sphere_conv_cut,isoval),...
	'FaceColor','blue',...
	'EdgeColor','none');
isonormals(Sphere_conv_cut,hisoxy)
toc

set(gcf,'Renderer','zbuffer'); lighting phong
set(hisoxy,'SpecularColorReflectance',0.2,'SpecularExponent',50)
lightangle(30,30);
view(24, 12)
axis off
grid off
axis equal

Proj = sum(Sphere_conv,3);




%%
figure('Color','white');
surf(gv,gv,Proj,'EdgeColor','none')
camlight left; lighting phong
colormap hot
axis off
grid off
view(-34, 58)

[X,Y] = meshgrid(gv);
LocError = exp(-(X.^2 + Y.^2)/Sigma^2);
figure('Color','white');
surf(gv,gv,LocError,'EdgeColor','none')
camlight left; lighting phong
colormap hot
axis off
grid off
view(-34, 58)

Virus_model = convn(Proj,LocError,'same');
figure('Color','white');
surf(gv,gv,Virus_model,'EdgeColor','none')
camlight left; lighting phong
colormap hot
axis off
grid off
view(-34, 58)





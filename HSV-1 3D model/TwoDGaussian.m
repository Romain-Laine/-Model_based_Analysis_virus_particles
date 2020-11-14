clc
clear all
close all


Step = 2; % nm


Sigma = 150; % nm

Max_xy = 4*Sigma;
gv = -Max_xy:Step:Max_xy;  % nm

%%

[X,Y] = meshgrid(gv);
LocError = exp(-(X.^2 + Y.^2)/(2*Sigma^2));
figure('Color','white');
surf(gv,gv,LocError,'EdgeColor','none');
camlight left; lighting phong

map = colormap(hot(512));
map(257:end,:) = [];
colormap(map);

axis off
grid off
view(-38, 38)



figure('Color','white');
surf(gv,gv,LocError);
camlight left; lighting phong

map = colormap(hot);
map(33:end,:) = [];
colormap(map);

axis off
grid off
view(-38, 38)



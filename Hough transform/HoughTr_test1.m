clear all 
close all
clc

N = 200;
r = 80; % nm
Phi = 2*pi*rand(N,1);

x0 = 100;
y0 = 20;

PixelSize = 10; % nm
Size = max(abs(x0),abs(y0)) + 2*r;

xy = [x0 + r.*cos(Phi) y0 + r.*sin(Phi)];

figure('Color','white')
plot(xy(:,1),xy(:,2),'+')
axis equal
grid on



x_grid = -Size:PixelSize:Size;
xy_round = PixelSize*round(xy/PixelSize);
Density_image = zeros(size(x_grid,2),size(x_grid,2));

%%
for i = 1:size(xy_round,1)
    Density_image((x_grid == xy_round(i,1)),(x_grid == xy_round(i,2))) = Density_image((x_grid == xy_round(i,1)),(x_grid == xy_round(i,2))) + 1;        
end

%%
h = figure('Color','white');
imshow(Density_image,[])

radius = round(r/PixelSize*[0.5, 1.5]);
[centers, radii] = imfindcircles(Density_image,radius,'ObjectPolarity', 'bright', 'Sensitivity',0.95);

disp(['Simulated values:',num2str(r),' ',num2str(y0),' ',num2str(x0)]);

calc_radius = radii*PixelSize;
calc_center = PixelSize*centers - PixelSize*(1+(size(x_grid,2)-1)/2);

disp(['Calculated values:',num2str(calc_radius),' ',num2str(calc_center)]);




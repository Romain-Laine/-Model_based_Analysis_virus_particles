function [ ImCircle, ImCentres, LogicalIm ] = CreateCircleImage( Image_size, Centres, Radii )

[x_grid, y_grid] = meshgrid(1:Image_size(2),1:Image_size(1));

n_circles = size(Radii,1);
ImCircle = zeros(Image_size(2),Image_size(1));
ImCentres = zeros(Image_size(2),Image_size(1));
LogicalIm = zeros(Image_size(2),Image_size(1));

disp('Creating circles...')
h_wait = waitbar(0,'Please wait while the circles are drawn...','name','Wait...');

for i = 1:n_circles
    ImCentres(round(Centres(i,2)),round(Centres(i,1)))= 1;
    LogicalIm = LogicalIm + ((x_grid - Centres(i,1)).^2 + (y_grid - Centres(i,2)).^2 < (1.1*Radii(i))^2);
    ImCircle = ImCircle + ((x_grid - Centres(i,1)).^2 + (y_grid - Centres(i,2)).^2 < (1.1*Radii(i))^2) - ((x_grid - Centres(i,1)).^2 + (y_grid - Centres(i,2)).^2 < (0.95*Radii(i))^2);
    waitbar(i/n_circles);
end

close(h_wait);
ImCircle = im2bw(ImCircle, 0);
ImCentres = im2bw(ImCentres, 0);
LogicalIm = im2bw(LogicalIm, 0);

end


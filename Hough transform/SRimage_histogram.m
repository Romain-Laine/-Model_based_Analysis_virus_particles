function [ h ] = SRimage_histogram(xy, PixelNumber, Size )
% SRimage_histogram displays localization as a density map through
% simple histogramming


xi = linspace(-1,1,PixelNumber)*Size;
xr = interp1(xi,1:numel(xi),xy(:,1),'nearest');
yr = interp1(xi,1:numel(xi),xy(:,2),'nearest');
LocEdge = find(isnan(xr) == 1 | isnan(yr) == 1);
xr(LocEdge) = [];
yr(LocEdge) = [];

Density_image = accumarray([yr xr],1);

h = figure('Color','white');
imshow(Density_image,[])


end


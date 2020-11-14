function SR_image = SRimage_histogram(xy, PixelNumber, PixelSize )
% SRimage_histogram displays localization as a density map through
% simple histogramming

SR_image = zeros(PixelNumber);

xi = (0:PixelNumber)*PixelSize;
xr = interp1(xi,1:numel(xi),xy(:,1),'nearest');
yr = interp1(xi,1:numel(xi),xy(:,2),'nearest');
LocEdge = find(isnan(xr) == 1 | isnan(yr) == 1);
xr(LocEdge) = [];
yr(LocEdge) = [];

hist_image = accumarray([yr xr],1);
SR_image(1:size(hist_image,1),1:size(hist_image,2)) = hist_image;

end


function [ ImCircle ] = CreateCircleImage( Image, Centre, Radius )

[x_grid, y_grid] = meshgrid(1:size(Image,2),1:size(Image,1));
ImCircle = ((x_grid - Centre(:,1)).^2 + (y_grid - Centre(:,2)).^2 < (1.1*Radius)^2) - ((x_grid - Centre(:,1)).^2 + (y_grid - Centre(:,2)).^2 < (0.9*Radius)^2);

% figure;
% imshow(ImCircle,[])

end


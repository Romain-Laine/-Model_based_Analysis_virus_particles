function ImPanel = DispListObject_Centroids(I, Centroids, PixelSize, ScaleBarSize)
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here
% Centroids are already converted in pixel (from nm)

CropSize = 250; % nm
CropSize = CropSize/PixelSize;
n_object = size(Centroids,1);

ImObject = cell(1, n_object);

for i=1:n_object
    xmin = max(round(Centroids(i,1) - CropSize),0);
    ymin = max(round(Centroids(i,2) - CropSize),0);
    rect = [xmin ymin 2*CropSize 2*CropSize];
    ImObject{i} = imcrop(I,rect);
    ImObject{i} = AddScaleBar(ImObject{i}, PixelSize, ScaleBarSize );
end


BorderSize = 5; % in pixels
PanelSize = BorderSize + ceil(sqrt(n_object))*(BorderSize + 2*CropSize);
ImPanel = ones(PanelSize, PanelSize, 3);
ImPanel = uint8(256*ImPanel);

ImNumbering = 1:(ceil(sqrt(n_object)))^2;
ImNumbering = reshape(ImNumbering, ceil(sqrt(n_object)), ceil(sqrt(n_object)));
ImNumbering(ImNumbering > n_object) = 0;


for i = 1:ceil(sqrt(n_object))       % line
    for j = 1:ceil(sqrt(n_object))   % column
        if ImNumbering(i,j) == 0
        else
        ImSize = size(ImObject{ImNumbering(i,j)});
        ImPanel((j-1)*(2*CropSize + BorderSize) + 1 + BorderSize:(ImSize(1)+(j-1)*(2*CropSize + BorderSize)) + BorderSize,...
            (i-1)*(2*CropSize + BorderSize) + 1 + BorderSize:(ImSize(2)+(i-1)*(2*CropSize + BorderSize) + BorderSize),1:3) = ImObject{ImNumbering(i,j)};
        end
    end
end

figure('Color','white','name','List of objects centered on loaded centroids');
imshow(ImPanel)
title(['Number of objects: ', num2str(n_object),' - Scale bar: ',num2str(ScaleBarSize),' nm'])

end


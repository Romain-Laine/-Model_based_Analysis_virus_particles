function ImPanel = DispListObject(I, ParamObject, ListObject, CropSize, PixelSize, ScaleBarSize)
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here

n_object = max(size(ListObject));
ImObject = cell(1, n_object);

for i=1:n_object
    xy0 = ParamObject(ListObject(i)).Centroid;
    xmin = max(round(xy0(1) - CropSize),0);
    ymin = max(round(xy0(2) - CropSize),0);
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

figure('Color','white','name','Particle panel');
imshow(ImPanel)
title(['Number of objects: ', num2str(n_object),' - Scale bar: ',num2str(ScaleBarSize),' nm'])




end


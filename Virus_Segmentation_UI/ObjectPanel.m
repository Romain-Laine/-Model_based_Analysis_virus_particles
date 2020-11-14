function ImPanel = ObjectPanel(Image, ParamObject, ListObject, PixelSize, MaxPanelSize, Scale_bar_size)
%UNTITLED10 Summary of this function goes here
%   Detailed explanation goes here

% MaxPanelSize = 100;  % In number of images 4 --> 4 x 4 image panels are created

Diameters = [ParamObject(ListObject).EquivDiameter];
CropSize = round(mean(Diameters));

if max(size(ListObject)) < MaxPanelSize^2
    ImPanel = DispListObject(Image, ParamObject, ListObject, CropSize, PixelSize, Scale_bar_size);
else
    while max(size(ListObject)) > MaxPanelSize^2
        ImPanel = DispListObject(Image, ParamObject, ListObject(1:MaxPanelSize^2),CropSize, PixelSize, Scale_bar_size);
        ListObject(1:MaxPanelSize^2) = [];
    end
    
    if isempty(ListObject) == 1
    else
        ImPanel = DispListObject(Image, ParamObject, ListObject, CropSize, PixelSize, Scale_bar_size);
    end
end


end


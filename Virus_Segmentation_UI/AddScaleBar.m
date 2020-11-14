function ImObject_Out = AddScaleBar( ImObject, PixelSize, ScaleBarSize )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

ScaleBarSize = round(ScaleBarSize/PixelSize); % in pixels

ImObject_Out = ImObject;
ImSize = size(ImObject_Out);

if (ImSize(1) > 1) && (ImSize(2) > ScaleBarSize - 1)
    ImObject_Out(ImSize(1)-1, (ImSize(2)-ScaleBarSize-1):ImSize(2)-1, 1:3) = 256;
end


end


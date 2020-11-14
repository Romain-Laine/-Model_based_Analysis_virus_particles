function rect = CropObject( ParamObject, n)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
Diameter = ParamObject(n).EquivDiameter;
Centroid = ParamObject(n).Centroid;
rect = [round(Centroid(1)-Diameter) round(Centroid(2)-Diameter) 2*Diameter 2*Diameter];


end


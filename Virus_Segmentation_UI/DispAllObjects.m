function Image = DispAllObjects(I, ParamObject, ListObject)
%UNTITLED8 Summary of this function goes here
%   Detailed explanation goes here

n_object = max(size(ListObject));
Image = zeros(size(I));

if isfield(ParamObject,'PixelList')
    for i=1:n_object
        Coord = ParamObject(ListObject(i)).PixelList;
        for j = 1:size(Coord,1)
            Image(Coord(j,2),Coord(j,1)) = 1;
        end
    end
else
    [x_grid, y_grid] = meshgrid(1:size(I,2),1:size(I,1));
    Centres = cat(1, ParamObject.Centroid);
    Radii = cat(1, ParamObject.EquivDiameter)/2;
    
    for i = 1:n_object
        Image = Image + ((x_grid - Centres(ListObject(i),1)).^2 + (y_grid - Centres(ListObject(i),2)).^2 < (1.1*Radii(ListObject(i)))^2);
    end
end


end
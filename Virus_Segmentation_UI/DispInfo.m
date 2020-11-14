function DispInfo( ParamObject, ObjectList )
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

r_names = {'Area','Eccentricity','EquivDiameter'};
c_name = {ObjectList};

PanelSize = ceil(sqrt(max(size(ObjectList))));

Info = struct2cell(ParamObject(ObjectList));
Info(5,:) = [];    % discard info about pixellist
Info(2,:) = [];    % discard info about centroid

ColSize = 75;
% TableWidth = 150 + ColSize*max(size(ObjectList));
TableWidth = 150 + ColSize*PanelSize;

figure('Color','white','Position',[100 100 TableWidth 100]);
uitable('Data',Info,'RowName', r_names,'ColumnName',c_name,...
    'Units','normalized','Position',[0 0 1 1]);



end


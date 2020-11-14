function [ X, Y, Area] = SRimage2Centroids( Image, PixelSize, Gamma, Otsu_level, Min_object_size )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

% Area is the list of areas of all objects

Dilate_size = 3; % pixels
Open_size = 5; % pixels


%% Read image

if size(Image,3) > 1
    Image(:,:,2:3) = [];
end

if isa(Image, 'uint8')
    n = 8;
    %     disp('File type: uint8');
elseif isa(Image, 'uint16')
    n = 16;
    %     disp('File type: uint16');
end

%% Apply some gamma
Image2 = imadjust(Image,[0 1],[0 1],Gamma);


%% Determine Otsu level and threhsold

level = graythresh(Image2);
% disp(['Otsu threshold: ', num2str(level*Otsu_level)]);
bw = im2bw(Image2,level*Otsu_level);
bw = bwareaopen(bw, Open_size); % Get rid of the small objects
se = strel('disk',Dilate_size);
bw2 = imdilate(bw,se);

% Get rid of the small bits
cc = bwconncomp(bw2); % find teh connected regions
stats = regionprops(cc, 'Area');
idx = find([stats.Area] > Min_object_size);
bw3 = ismember(labelmatrix(cc), idx);
BW_outer = bwmorph(bw3,'remove');


Image_display = 2^n*double(Image2)/max(double(Image2(:)));
if isa(Image, 'uint8')
    Image_display = uint8(Image_display);
elseif isa(Image, 'uint16')
    Image_display = uint16(Image_display);
end

ColorImage = cat(3,Image_display, zeros(size(Image)), (2^n-1)*BW_outer);
figure('Color','white','name','Object segmentation image');
imshow(ColorImage);

%%

% Calculate area
Area = cell2mat(struct2cell(regionprops(bw3, 'Area'))); % compute the are of each object
Area = Area*PixelSize^2; % Area is calculated in nm^2 !!

figure('Color','white','name','Area histogram');
histogram(Area)
xlabel 'Area (nm^2)'
ylabel 'Occurences'
title(['Number of centroids found: ', num2str(length(Area))])

Centroid = cell2mat(reshape(struct2cell(regionprops(bw3, 'Centroid')),size(Area,2),size(Area,1)));  % in pixels
Centroid = PixelSize*Centroid; % in nm

if isempty(Centroid)
    X = zeros(0,0);
    Y = zeros(0,0);
    disp('WARNING: No centroids were found !');

else
    X = Centroid(:,1);
    Y = Centroid(:,2);
    disp(['Number of centroids found: ', num2str(size(Centroid,1))]);

end


end


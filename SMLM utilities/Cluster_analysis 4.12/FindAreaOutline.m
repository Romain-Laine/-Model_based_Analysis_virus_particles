function [ Mask ] = FindAreaOutline(PointImage, Analysis_window, SR_image_pixelSize )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here


Otsu_threshold = 0.1; % works well for the point images generated in the case of cluster analysis
Close_radius_size = 50;
Dilate_radius = ceil(2*Analysis_window/SR_image_pixelSize); % of the order of the analysis 2*analysis window 

Border = 2*(Close_radius_size + Dilate_radius);

bw = im2bw(PointImage,Otsu_threshold);
large_bw = false(size(PointImage,1)+2*Border, size(PointImage,2) + 2*Border);
large_bw((Border+1):(size(PointImage,1) + Border), (Border+1):(size(PointImage,2) + Border)) = bw;
% 
% figure;
% imshow(bw)
% figure;
% imshow(large_bw)


%% Dilate
se = strel('disk',Dilate_radius);
large_bw = imdilate(large_bw,se);

%% Close the image
se = strel('disk',Close_radius_size);
large_bw = imclose(large_bw, se);
% figure;
% imshow(large_bw)

Mask = large_bw((Border+1):(size(PointImage,1) + Border), (Border+1):(size(PointImage,2) + Border));
% figure;
% imshow(Mask)



end


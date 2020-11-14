function [ Image_out ] = AddCentroid2Image( Centroids, Image_in )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

n_row = size(Image_in,1);
n_col = size(Image_in,2);


Blank_image = zeros(n_row,n_col);
CentroidImage = Blank_image;

Centroids = floor(Centroids);

for i = 1:size(Centroids,1)
    CentroidImage(Centroids(i,2),Centroids(i,1)) = 1;
end


n_dilate = 1;
se = strel('disk',max(1,round(abs(n_dilate))));
CentroidImage = imdilate(CentroidImage, se);

Image_out = Image_in + uint8(256*cat(3,Blank_image,CentroidImage, Blank_image));

end


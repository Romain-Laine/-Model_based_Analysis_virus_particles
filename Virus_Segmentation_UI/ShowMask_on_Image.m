function [ Image_out ] = ShowMask_on_Image( Mask,  Image_in )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

numrows = size(Image_in,1);
numcols = size(Image_in,2);

Blank_image = zeros(numrows,numcols);

Mask = bwmorph(Mask,'remove');
Image_out = Image_in + uint8(256*cat(3,Blank_image,Blank_image, Mask));


end


function [ ColorImage ] = Create2Cimage( Image_RC, Image_GC, n )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

Image_height = size(Image_RC,1);
Image_width = size(Image_RC,2);

Image_RC = uint8((2^n-1)*Image_RC/max(Image_RC(:)));
Image_GC = uint8((2^n-1)*Image_GC/max(Image_GC(:)));

ColorImage = uint8(zeros(Image_height, Image_width, 3));
ColorImage(:,:,1) = Image_RC;
ColorImage(:,:,2) = Image_GC;



end


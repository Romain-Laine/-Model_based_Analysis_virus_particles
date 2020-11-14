function [ ColorImage ] = Grey2Color( GreyImage )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here


n = 8;
Cmap = colormap(hot(2^n));
Cmap = cat(1, [0 0 0],Cmap);
Cmap(end,:) = [];

ColorImage = uint8((2^n-1)*GreyImage);
ColorImage = ind2rgb(ColorImage, Cmap);
ColorImage = uint8((2^n-1)*ColorImage);

end


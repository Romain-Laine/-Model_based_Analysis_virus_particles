function HistParam(ParamObject)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

Sizes = [ParamObject.Eccentricity];

figure('Color','white');
hist(Sizes,100);
xlabel 'Eccentricity'
ylabel 'Occurences'
title 'Histogram of eccentricity'
xlim([0.01 1])

end


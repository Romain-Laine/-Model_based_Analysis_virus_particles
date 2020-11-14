function HistSize( ParamObject )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

Sizes = [ParamObject.Area];

figure('Color','white');
subplot(2,1,1)
hist(Sizes,100);
xlabel 'Number of pixel'
ylabel 'Occurences'
title 'Histogram of particle sizes (in pixel)'
subplot(2,1,2)
hist(Sizes,100);
xlabel 'Number of pixel'
ylabel 'Occurences'
ylim([0 10])

end


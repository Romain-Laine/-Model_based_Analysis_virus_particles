function  Mask = OtsuThreshold(I, Sensitivity)

smoothing = 5;
closing = 5;

se_smoothing = strel('disk',max(1,round(abs(smoothing))));
se_closing = strel('disk',max(1,round(abs(closing))));

hsize = 5;
sigma = 1;
filter = fspecial('gaussian', hsize, sigma);

I = imfilter(I, filter,'replicate','same');
% I = medfilt2(I,[3 3]); % remove some noise

I = imadjust(I);

Level = graythresh(I);
Level = max(Level/Sensitivity,0);
Level = min(Level/Sensitivity,1);
Mask = im2bw(I,Level);
Mask = bwmorph(Mask,'clean');

% Closes
n_close = 1;
for i = 1:n_close
    Mask = imclose(Mask,se_closing);
end

% Smoothes
Mask = imerode(Mask,se_smoothing);
Mask = imdilate(Mask,se_smoothing);


end
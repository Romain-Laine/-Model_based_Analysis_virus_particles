function [ hFigure ] = COM_display( ParticleInfo , ScaleBarSize, Gamma, I_cutoff)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

Xc = ParticleInfo(:,1);
Yc = ParticleInfo(:,2);
Framec = ParticleInfo(:,3);
ADCc = ParticleInfo(:,4);
nObject = ParticleInfo(:,5);

R = sqrt(Xc.^2+Yc.^2);
r = mean(R);
t = std(R);   % the standard deviation is half the thickness of the layer

hFigure = figure('Color','white','name','Centre of mass analysis');
subplot(2,6,1:3)
plot(Xc,Yc,'.','MarkerSize',15)
hold on
plot(0,0,'r+','MarkerSize',10,'LineWidth',2)
xlabel 'nm'
ylabel 'nm'
axis equal
grid on

ang = 0:0.01:2*pi;
xp = (r+t)*cos(ang);
yp = (r+t)*sin(ang);
plot(xp,yp,'r','LineWidth',2);

xp = (r-t)*cos(ang);
yp = (r-t)*sin(ang);
plot(xp,yp,'r','LineWidth',2);

formatSpec = '%10.1f';
title(['Radius = ',num2str(r,formatSpec),' nm & thickness = ',num2str(2*t,formatSpec),' nm'])

% Pixel size of the final super-resolution image
PixelSize = 10; % nm
Min_display_X = (ceil(abs(min(Xc))/PixelSize))*PixelSize;
Min_display_Y = (ceil(abs(min(Yc))/PixelSize))*PixelSize;

LocInfo = cat(2,Xc + Min_display_X, Yc + Min_display_Y, Framec, ADCc);
GreyImage = Loc2SRimage(LocInfo, PixelSize);
GreyImage = imadjust(GreyImage,[0 I_cutoff],[0 1],Gamma);
ColorImage = Grey2Color(GreyImage);

% Write the image in the base workspace before the scale bar is added
% assignin('base', 'ColorImage', ColorImage);
% Filename = 'C:\Users\rfl30\HSV-1 virus particle\Paper drafts\Images\Single particles\Aligned particles\AlignedParticle.png';
% imwrite(ColorImage,Filename,'png');


ColorImage = AddScaleBar( ColorImage, PixelSize, ScaleBarSize);

subplot(2,6,7:9)
imshow(ColorImage)
title(['Scale bar: ',num2str(ScaleBarSize),' nm (',num2str(PixelSize),' nm / pixel)'])

subplot(2,6,[5 6 11 12])
hist(R,round(0.5*sqrt(max(size(R)))))
xlabel 'radius (nm)'
ylabel 'Number of localizations'
title([num2str(max(size(R))),' localizations from ', num2str(max(nObject)),' particles'])



end


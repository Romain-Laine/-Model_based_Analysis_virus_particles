% circ2dGaussVirusEstParams
%
% V1-1, EJR, November 2013
%
% PURPOSE:
%    Estimates observable RDF parameters of virus localisation data
%
% METHOD
%    Takes Virus and Localisation Microscopy parameters as input
%    Input Parameters are: radVP  (virus protein radius) 
%                    sigPos (rms distance of localisation from VP position)
%    Uses a grid (1 nm squares - fine enough?) to determine observables
%    Outputs are the observable parameters of the localisation RDF: 
%     Mean and Std of radial distance of localisations. 
%
% COMMENTS:
%    I also have a 1-D integral which might be easier to integrate
%    This script is more of a physical simulation


function [estRadMean, estRadStd] = circ2dGaussVirusGetParams(inputs)

radVP = inputs(1);
sigPos = inputs(2);

% Find observed parameters, to enable error determination. 
% OR import as an Anonymous Function
obsMeanRad = mean( evalin('base','radLocs') );
obsStdRad  =  std( evalin('base','radLocs') );


% Create a suitably large, fine simulation grid
gridRadius = ceil ( obsMeanRad + obsStdRad*3 ); 
gridSpacing = 1; % nm - initial guess; could improve

blurGridDiam = 6*ceil(obsStdRad);
gridBlur = fspecial('gaussian', blurGridDiam, abs(sigPos) ); 
% gridBlur = evalin('base','gridBlur');

% Create an underlying structure to blur:
im = zeros( (2*gridRadius/gridSpacing + 1) );
radius = radVP;
xc = gridRadius +1;
yc = gridRadius +1;
value = 1;

% -------- Midpoint circle algorithm
% Draw a circle in a matrix using the integer midpoint circle algorithm
% Does not miss or repeat pixels
% Created by : Peter Bone
% Created : 19th March 2007
% % im = MidpointCircle(im, radius, xc, yc, value)

xc = int16(xc);
yc = int16(yc);

x = int16(0);
y = int16(radius);
d = int16(1 - radius);

im(xc, yc+y) = value;
im(xc, yc-y) = value;
im(xc+y, yc) = value;
im(xc-y, yc) = value;

while ( x < y - 1 )
    x = x + 1;
    if ( d < 0 ) 
        d = d + x + x + 1;
    else 
        y = y - 1;
        a = x - y + 1;
        d = d + a + a;
    end
    im( x+xc,  y+yc) = value;
    im( y+xc,  x+yc) = value;
    im( y+xc, -x+yc) = value;
    im( x+xc, -y+yc) = value;
    im(-x+xc, -y+yc) = value;
    im(-y+xc, -x+yc) = value;
    im(-y+xc,  x+yc) = value;
    im(-x+xc,  y+yc) = value;
end
 
% ----------- End of midpoint circle algorithm

imBlurred = conv2(im, gridBlur,'same');
sumImBlur = sum(imBlurred(:));

[XX,YY] = meshgrid(-gridRadius:gridRadius,-gridRadius:gridRadius);

myRads = sqrt((XX.^2 + YY.^2));

estRadMean = sum(sum( (myRads.*imBlurred) )) / sumImBlur;
estRadStd  = sqrt( sum(sum((myRads-estRadMean).^2.*imBlurred))/sumImBlur);

assignin('base','estRadMean',estRadMean);
assignin('base','estRadStd',estRadStd);


estRadErr = (estRadMean - obsMeanRad)*1;
estSigErr = (estRadStd  - obsStdRad)*1;

end
function [ SR_image ] = SRdisplay_rapidSTORM( xy, PixelSize, SR_imageSize )
%SRDISPLAY_RAPIDSTORM Summary of this function goes here
%   Detailed explanation goes here

% Transform coordinates in unit of PixelSize
xy = xy/PixelSize;

% Number of localizations
n_loc = size(xy,1);
SR_image = zeros(SR_imageSize);

for i = 1:n_loc
    x0 = round(xy(i,2));   % 0-based
    y0 = round(xy(i,1));   % 0-based
    SR_image(x0+1,y0+1) = SR_image(x0+1,y0+1) + W_rS(xy(i,2) - x0    )*W_rS(xy(i,1) - y0);       % (x0,y0)
    SR_image(x0+2,y0+1) = SR_image(x0+2,y0+1) + W_rS(xy(i,2) - (x0+1))*W_rS(xy(i,1) - y0);       % (x0+1,y0)
    SR_image(x0+1,y0+2) = SR_image(x0+1,y0+2) + W_rS(xy(i,2) - x0    )*W_rS(xy(i,1) - (y0+1));   % (x0,y0+1)
    SR_image(x0+2,y0+2) = SR_image(x0+2,y0+2) + W_rS(xy(i,2) - (x0+1))*W_rS(xy(i,1) - (y0+1));   % (x0+1,y0+1)
    
    if x0 > 0
        SR_image(x0  ,y0+1) = SR_image(x0  ,y0+1) + W_rS(xy(i,2) - (x0-1))*W_rS(xy(i,1) - y0);       % (x0-1,y0)
        SR_image(x0  ,y0+2) = SR_image(x0  ,y0+2) + W_rS(xy(i,2) - (x0-1))*W_rS(xy(i,1) - (y0+1));   % (x0-1,y0+1)
    end
    
    if y0 > 0
        SR_image(x0+1,y0  ) = SR_image(x0+1,y0  ) + W_rS(xy(i,2) - x0    )*W_rS(xy(i,1) - (y0-1));   % (x0,y0-1)
        SR_image(x0+2,y0  ) = SR_image(x0+2,y0  ) + W_rS(xy(i,2) - (x0+1))*W_rS(xy(i,1) - (y0-1));   % (x0+1,y0-1)
    end
    
    if (x0>0) && (y0>0)
        SR_image(x0  ,y0  ) = SR_image(x0  ,y0  ) + W_rS(xy(i,2) - (x0-1))*W_rS(xy(i,1) - (y0-1));   % (x0-1,y0-1)
    end
    
end


end


function W = W_rS(x)

if abs(x)>1
    W = 0;
else
    W = 1-abs(x);
end

end

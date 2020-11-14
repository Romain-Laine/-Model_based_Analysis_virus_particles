function [ SR_Image ] = Loc2SRimage( LocInfo, PixelSize)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

Method = 2;
% 0 = histogramming
% 1 = histogramming with some nearest neighbour smoothing
% 2 = histogramming with pixelization identical to rapidSTORM

X = LocInfo(:,1); % in nm
Y = LocInfo(:,2); % in nm
% Frame = LocInfo(:,3);
ADC = LocInfo(:,4);

PixelNumber = ceil(max(max(X),max(Y))/PixelSize); % number of pixels in the SR image
disp('Creating the STORM image from the localization file');

tic
if Method == 0
    xi = (0:PixelSize:(PixelNumber - 1)*PixelSize);
    xr = interp1(xi,1:numel(xi),X,'nearest');
    yr = interp1(xi,1:numel(xi),Y,'nearest');
    LocEdge = find(isnan(xr) == 1 | isnan(yr) == 1);
    xr(LocEdge) = [];
    yr(LocEdge) = [];
    ADC(LocEdge) = [];
    
    SR_Image = accumarray([yr xr],ADC,[PixelNumber PixelNumber]);
    
elseif Method == 1
    X = X/PixelSize;
    Y = Y/PixelSize;
    n_loc = size(X,1);
    SR_Image = zeros(PixelNumber + 1); % an extra ipxel is necessary to take into account that pixel pixel (1,1) is actually (0,0)nm
    
    for i = 1:n_loc
        
        % The max(..., 0) is a failsafe for when X or Y is negative
        x1 = max(floor(X(i)),0);
        y1 = max(floor(Y(i)),0);
        x2 = max(ceil(X(i)),0);
        y2 = max(ceil(Y(i)),0);
        
        % Calculate the distance with respect to the centre of all adjacent pixels
        D11 = sqrt((x1-X(i))^2 + (y1-Y(i))^2);
        D12 = sqrt((x2-X(i))^2 + (y1-Y(i))^2);
        D21 = sqrt((x1-X(i))^2 + (y2-Y(i))^2);
        D22 = sqrt((x2-X(i))^2 + (y2-Y(i))^2);
        
        % Calculate the distance with respect to the centre of all adjacent pixels
        %         D11 = max(0,1 - D11);
        %         D12 = max(0,1 - D12);
        %         D21 = max(0,1 - D21);
        %         D22 = max(0,1 - D22);
        
        % Transform the distances into a gaussian to favour the pixel in the centre
        G = 2;
        D11 = exp(-G*D11^4);
        D12 = exp(-G*D12^4);
        D21 = exp(-G*D21^4);
        D22 = exp(-G*D22^4);
        
        S = D11 + D12 + D21 + D22;
        
        if S == 0
            SR_Image(y1+1,x1+1) = ADC(i); % case where the localization is integer values like (1,1) or (4,4), quite unlikely
        else
            % The proportion of signal is inversely proportional to the distance to the centre of the pixel
            SR_Image(y1+1,x1+1) = SR_Image(y1+1,x1+1) + ADC(i)*D11/S;
            SR_Image(y2+1,x1+1) = SR_Image(y2+1,x1+1) + ADC(i)*D12/S;
            SR_Image(y1+1,x2+1) = SR_Image(y1+1,x2+1) + ADC(i)*D21/S;
            SR_Image(y2+1,x2+1) = SR_Image(y2+1,x2+1) + ADC(i)*D22/S;
        end
        
    end
    
    
elseif Method == 2
    border = 2;   % pixels
    % The max is the max of the radius
%     xy_max = max(sqrt(X.^2 + Y.^2))
    X(:,1) = X(:,1) + PixelSize*border;
    Y(:,1) = Y(:,1) + PixelSize*border;
    
    % Calculate the image size with some pixels on the border
    SR_imageSize = ceil(max(max(X),max(Y))/PixelSize) + border;
    SR_Image = SRdisplay_rapidSTORM( cat(2,X,Y), PixelSize, SR_imageSize);
    
end
SR_Image = SR_Image/max(max(SR_Image)); % normalize
toc

end


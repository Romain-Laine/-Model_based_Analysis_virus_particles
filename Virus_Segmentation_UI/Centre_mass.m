function ParticleInfo = Centre_mass(LocInfo, Centroids, PixelSize, Radius_cutoff)

X = LocInfo (:,1); % in nm
Y = LocInfo (:,2); % in nm
Frame = LocInfo(:,3);
ADC = LocInfo(:,4);
n_object = size(Centroids,1);

C_nm = (Centroids - 1)*PixelSize;  % Correcting for the pixel offset of the centroid with respect to the SR image

Xc = 0;
Yc = 0;
ADCc = 0;
Framec = 0;
nObject = 0;


for i = 1:n_object
    R = sqrt((X-C_nm(i,1)).^2 + (Y-C_nm(i,2)).^2);
    Loc = find(R <= Radius_cutoff);
    Xc = cat(1,Xc,X(Loc)-C_nm(i,1));
    Yc = cat(1,Yc,C_nm(i,2)-Y(Loc));
    Framec = cat(1,Framec,Frame(Loc));
    ADCc = cat(1,ADCc,ADC(Loc));
    nObject = cat(1,nObject,i*ones(size(Loc,1),1));
end

Xc(1) = [];
Yc(1) = [];
ADCc(1) = [];
Framec(1) = [];
nObject(1) = [];

ParticleInfo = cat(2,Xc,Yc,Framec,ADCc,nObject);







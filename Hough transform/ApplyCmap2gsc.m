function [ RGB_Im ] = ApplyCmap2gsc( GSC_Image, Cmap )

RGB_Im = uint8(zeros(size(GSC_Image,1),size(GSC_Image,2),3));
GSC_Image = uint8(255*GSC_Image/max(max(GSC_Image)));

for i = 1:size(GSC_Image,1)
    for j = 1:size(GSC_Image,2)
        RGB_Im(i,j,1) = 255*Cmap(GSC_Image(i,j)+1,1);
        RGB_Im(i,j,2) = 255*Cmap(GSC_Image(i,j)+1,2);
        RGB_Im(i,j,3) = 255*Cmap(GSC_Image(i,j)+1,3);
    end
end




end


function [ imgOut ] = pivMinMaxFilter( imgIn,mask,kernel,mindiff )
% PIVMINMAXFILTER Applies min-max filtering to an image (

% create image with NaN borders and with masked pixels set to NaN's
imgIn = single(imgIn);
imgIn(mask==0) = NaN;
auxImg = zeros(size(imgIn)+[1 1]*(kernel-1))+NaN;
auxImg((kernel-1)/2:end-(kernel+1)/2,(kernel-1)/2:end-(kernel+1)/2) = imgIn;
auxMin = 0*auxImg;
auxMax = 0*auxImg;
auxMin2 = 0*auxImg;
auxMax2 = 0*auxImg;

% get minima and maxima
for ky = (kernel-1)/2:size(auxImg,1)-(kernel+1)/2
    for kx = (kernel-1)/2:size(auxImg,2)-(kernel+1)/2
        aux = auxImg(ky-(kernel-1)/2+1:ky+(kernel-1)/2+1,kx-(kernel-1)/2+1:kx+(kernel-1)/2)+1;
        aux = reshape(aux,numel(aux),1);
        aux = aux(~isnan(aux));
        if numel(aux)>1
            auxMin(ky,kx) = min(aux);
            auxMax(ky,kx) = max(aux);
        else
            auxMin(ky,kx) = NaN;
            auxMax(ky,kx) = NaN;
        end
    end
end

% filter minima and maxima
for ky = (kernel-1)/2:size(auxImg,1)-(kernel+1)/2
    for kx = (kernel-1)/2:size(auxImg,2)-(kernel+1)/2
        aux = auxMin(ky-(kernel-1)/2+1:ky+(kernel-1)/2+1,kx-(kernel-1)/2+1:kx+(kernel-1)/2)+1;
        aux = reshape(aux,numel(aux),1);
        aux = aux(~isnan(aux));
        if numel(aux)>1
            auxMin2(ky,kx) = sum(aux,1)/numel(aux);
        else
            auxMin2(ky,kx) = NaN;
        end
        aux = auxMax(ky-(kernel-1)/2+1:ky+(kernel-1)/2+1,kx-(kernel-1)/2+1:kx+(kernel-1)/2)+1;
        aux = reshape(aux,numel(aux),1);
        aux = aux(~isnan(aux));
        if numel(aux)>1
            auxMax2(ky,kx) = sum(aux,1)/numel(aux);
        else
            auxMax2(ky,kx) = NaN;
        end
    end
end

diff = auxMax2 - auxMin2;
diff(diff < mindiff) = mindiff;
auxMax2 = auxMin2 + diff;

auxMin2 = auxMin2((kernel-1)/2:end-(kernel+1)/2,(kernel-1)/2:end-(kernel+1)/2);
auxMax2 = auxMax2((kernel-1)/2:end-(kernel+1)/2,(kernel-1)/2:end-(kernel+1)/2);

imgOut = (imgIn - auxMin2)./(auxMax2-auxMin2) * 128;
end


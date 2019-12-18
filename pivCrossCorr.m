function [pivData,ccPeakIm] = pivCrossCorr(exIm1,exIm2,pivData,pivPar)
% pivCrossCorr - cross-correlates two images to find object displacement between them
%
% Usage:
% [pivData,ccPeakIm] = pivCrossCorr(exIm1,exIm2,pivData,pivPar)
%
% Inputs:
%    Exim1,Exim2 ... pair of expanded images (use pivInterrogate for their generation)
%    pivData ... (struct) structure containing more detailed results. Following fields are required (use
%             pivInterrogate for generating them):
%        Status ... matrix describing status of velocity vectors (for values, see Outputs section)
%        iaX, iaY ... matrices with centers of interrogation areas
%        iaU0, iaV0 ... mean shift of IA's
%    pivPar ... (struct) parameters defining the evaluation. Following fields are considered:
%       ccRemoveIAMean ... if =0, do not remove IA's mean before cross-correlation; if =1, remove the mean; if
%                          in between, remove the mean partially
%       ccMaxDisplacement ... maximum allowed displacement to accept cross-correlation peak. This parameter is 
%                             a multiplier; e.g. setting ccMaxDisplacement = 0.6 means that the cross-
%                             correlation peak must be located within [-0.6*iaSizeX...0.6*iaSizeX,
%                             -0.6*iaSizeY...0.6*iaSizeY] from the zero displacement. Note: IA offset is not 
%                             included in the displacement, hence real displacement can be larger if 
%                             ccIAmethod is any other than 'basic'.
%       ccWindow ... windowing function p. 389 in ref [1]). Possible values are 'uniform', 'gauss',
%             'parzen', 'hanning', 'welch'
%       ccCorrectWindowBias ... Set whether cc is corrected for bias due to IA windowing.
%            - default value: true ( will correct peak )
%       ccMethod ... set methods for finding cross-correlation of interrogation areas. Possible values are
%              'fft' (use Fast Fourrier Transform' and 'dcn' (use discrete convolution). Option 'fft' is more
%              suitable for initial guess of velocity. Option 'dcn' is suitable for final iterations of the
%              velocity field, if the displacement corrections are already small.
%       ccMaxDCNdist ... Defines maximum displacement, for which the correlation is computed by DCN method 
%              (apllies only for ccMethod = 'dcn'). 
% Outputs:
%    pivData  ... (struct) structure containing more detailed results. Following fields are added or updated:
%        X, Y, U, V ... contains velocity field
%        Status ... matrix with statuis of velocity vectors (uint8). Bits have this coding:
%            1 ... masked (set by pivInterrogate)
%            2 ... cross-correlation failed (set by pivCrossCorr)
%            4 ... peak detection failed (set by pivCrossCorr)
%            8 ... indicated as spurious by median test (set by pivValidate)
%           16 ... interpolated (set by pivReplaced)
%           32 ... smoothed (set by pivSmooth)
%        ccPeak ... table with values of cross-correlation peak
%        ccPeakSecondary ... table with values of secondary cross-correlation peak (maximum of
%            crosscorrelation, if 5x5 neighborhood of primary peak is removed)
%        ccStd1, ccStd2 ... tables with standard deviation of pixel values in interrogation area, for the
%            first and second image in the image pair
%        ccMean1, ccMean2 ... tables with mean of pixel values in interrogation area, for the
%            first and second image in the image pair
%        ccFailedN ... number of vectors for which cross-correlation failed
%            at distance larger than ccMaxDisplacement*(iaSizeX,iaSizeY) )
%        ccSubpxFailedN ... number of vectors for which subpixel interpolation failed
%      - Note: fields iaU0 and iaV0 are removed from pivData
%    ccPeakIm ... expanded image containing cross-correlation functions for all IAs (normalized by .ccPeak)
% 
% Important local variables:
%    failFlag ... contains value of status elements of the vector being processed
%    Upx, Vpx ... rough position of cross-correlation peak (before subpixel interpolation, in integer 
%                     number of pixels)
%
%%
% This subroutine is a part of
%
% =========================================
%               PIVsuite
% =========================================
%
% PIVsuite is a set of subroutines intended for processing of data acquired with PIV (particle image
% velocimetry) within Matlab environment.
%
% Written by Jiri Vejrazka, Institute of Chemical Process Fundamentals, Prague, Czech Republic
%
% For the use, see files example_XX_xxxxxx.m, which acompany this file. PIVsuite was tested with
% Matlab 8.2 (R2013b).
%
% In the case of a bug, please, contact me: vejrazka (at) icpf (dot) cas (dot) cz
%
%
% Requirements:
%     Image Processing Toolbox
%         (required only if pivPar.smMethod is set to 'gaussian')
%
%     inpaint_nans.m
%         subroutine by John D'Errico, available at http://www.mathworks.com/matlabcentral/fileexchange/4551
%
%     smoothn.m
%         subroutine by Damien Garcia, available at
%         http://www.mathworks.com/matlabcentral/fileexchange/274-smooth
%
% Credits:
%    PIVsuite is a redesigned version of a part of PIVlab software [3], developped by W. Thielicke and
%    E. J. Stamhuis. Some parts of this code are copied or adapted from it (especially from its
%    piv_FFTmulti.m subroutine).
%
%    PIVsuite uses 3rd party software:
%        inpaint_nans.m, by J. D'Errico, [2]
%        smoothn.m, by Damien Garcia, [5]
%
% References:
%   [1] Adrian & Whesterweel, Particle Image Velocimetry, Cambridge University Press, 2011
%   [2] John D'Errico, inpaint_nans subroutine, http://www.mathworks.com/matlabcentral/fileexchange/4551
%   [3] W. Thielicke and E. J. Stamhuid, PIVlab 1.31, http://pivlab.blogspot.com
%   [4] Raffel, Willert, Wereley & Kompenhans, Particle Image Velocimetry: A Practical Guide. 2nd edition,
%       Springer, 2007
%   [5] Damien Garcia, smoothn subroutine, http://www.mathworks.com/matlabcentral/fileexchange/274-smooth
%%


%% 0. Initialization

U = pivData.U;
V = pivData.V;
status = pivData.Status;
ccPeak = U;                  % will contain peak levels
ccPeakSecondary = U;         % will contain level of secondary peaks
iaSizeX = pivPar.iaSizeX;
iaSizeY = pivPar.iaSizeY;
iaNX = size(pivData.X,2);
iaNY = size(pivData.X,1);
ccStd1 = pivData.U+NaN;
ccStd2 = pivData.U+NaN;
ccMean1 = pivData.U+NaN;
ccMean2 = pivData.U+NaN;


% initialize "expanded image" for storing cross-correlations
ccPeakIm = exIm1 + NaN;   % same size as expanded images

% peak position is shifted by 1 or 0.5 px, depending on IA size
if rem(iaSizeX,2) == 0
    ccPxShiftX = 1;
else
    ccPxShiftX = 0.5;
end
if rem(iaSizeY,2) == 0
    ccPxShiftY = 1;
else
    ccPxShiftY = 0.5;
end

%% 1. Create windowing function W and loss-of-correlation function F
% (ref. [1], Table 8.1, p. 390)
auxX = ones(iaSizeY,1)*(-(iaSizeX-1)/2:(iaSizeX-1)/2);
auxY = (-(iaSizeY-1)/2:(iaSizeY-1)/2)'*ones(1,iaSizeX);
EtaX = auxX/iaSizeX;
EtaY = auxY/iaSizeY;
KsiX = 2*EtaX;
KsiY = 2*EtaY;
if ~isfield(pivPar,'ccWindow') || strcmpi(pivPar.ccWindow,'uniform')
    W = ones(iaSizeY,iaSizeX);
    F = (1-abs(KsiX)).*(1-abs(KsiY));
elseif strcmpi(pivPar.ccWindow,'parzen')
    % window function W
    W = (1-2*abs(EtaX)).*(1-2*abs(EtaY));
    % loss of correlation function F
    auxFx = auxX + NaN; % initialization
    auxFy = auxX + NaN;
    auxOK = logical(abs(KsiX)<=1/2);
    auxFx(auxOK) = 1-6*KsiX(auxOK).^2+6*abs(KsiX(auxOK)).^3;
    auxFx(~auxOK) = 2-6*abs(KsiX(~auxOK))+6*KsiX(~auxOK).^2-2*abs(KsiX(~auxOK)).^3;
    auxOK = logical(abs(KsiY)<=1/2);
    auxFy(auxOK) = 1-6*KsiY(auxOK).^2+6*abs(KsiY(auxOK)).^3;
    auxFy(~auxOK) = 2-6*abs(KsiY(~auxOK))+6*KsiY(~auxOK).^2-2*abs(KsiY(~auxOK)).^3;
    F = auxFx.*auxFy;
elseif strcmpi(pivPar.ccWindow,'Hanning')
    W = (1/2+1/2*cos(2*pi*EtaX)).*(1/2+1/2*cos(2*pi*EtaY));
    F = (2/3*(1-abs(KsiX)).*(1+1/2*cos(2*pi*KsiX))+1/2/pi*sin(2*pi*abs(KsiX))) .*...
        (2/3*(1-abs(KsiY)).*(1+1/2*cos(2*pi*KsiY))+1/2/pi*sin(2*pi*abs(KsiY)));
elseif strcmpi(pivPar.ccWindow,'Welch')
    W = (1-(2*EtaX).^2).*(1-(2*EtaY).^2);
    F = (1-5*KsiX.^2+5*abs(KsiX).^3-abs(KsiX).^5) .* ...
        (1-5*KsiY.^2+5*abs(KsiY).^3-abs(KsiY).^5);
elseif strcmpi(pivPar.ccWindow,'Gauss')
    W = exp(-8*EtaX.^2).*exp(-8*EtaY.^2);
    F = exp(-4*KsiX.^2).*exp(-4*KsiY.^2);
elseif strcmpi(pivPar.ccWindow,'Gauss1')
    W = exp(-8*(EtaX.^2+EtaY.^2)) - exp(-2);
    W(W<0) = 0;
    W = W/max(max(W));
    F = NaN;
elseif strcmpi(pivPar.ccWindow,'Gauss2')
    W = exp(-16*(EtaX.^2+EtaY.^2)) - exp(-4);
    W(W<0) = 0;
    W = W/max(max(W));
    F = NaN;
elseif strcmpi(pivPar.ccWindow,'Gauss0.5')
    W = exp(-4*(EtaX.^2+EtaY.^2)) - exp(-1);
    W(W<0) = 0;
    W = W/max(max(W));
    F = NaN;
elseif strcmpi(pivPar.ccWindow,'Nogueira')
    W = 9*(1-4*abs(EtaX)+4*EtaX.^2).*(1-4*abs(EtaY)+4*EtaY.^2);
    F = NaN;
elseif strcmpi(pivPar.ccWindow,'Hanning2')
    W = (1/2+1/2*cos(2*pi*EtaX)).^2.*(1/2+1/2*cos(2*pi*EtaY)).^2;
    F = NaN;
elseif strcmpi(pivPar.ccWindow,'Hanning4')
    W = (1/2+1/2*cos(2*pi*EtaX)).^4.*(1/2+1/2*cos(2*pi*EtaY)).^4;
    F = NaN;
end
% Limit F to not be too small
F(F<0.5) = 0.5;

%% 2. Cross-correlate expanded images and do subpixel interpolation
% loop over interrogation areas
for kx = 1:iaNX
    for ky = 1:iaNY
        failFlag = status(ky,kx);
        % if not masked, get individual interrogation areas from the expanded images
        if failFlag == 0
            imIA1 = exIm1(1+(ky-1)*iaSizeY:ky*iaSizeY,1+(kx-1)*iaSizeX:kx*iaSizeX);
            imIA2 = exIm2(1+(ky-1)*iaSizeY:ky*iaSizeY,1+(kx-1)*iaSizeX:kx*iaSizeX);
            % remove IA mean
            auxMean1 = mean2(imIA1);
            auxMean2 = mean2(imIA2);
            imIA1 = imIA1 - pivPar.ccRemoveIAMean*auxMean1;
            imIA2 = imIA2 - pivPar.ccRemoveIAMean*auxMean2;
            % apply windowing function
            imIA1 = imIA1.*W;
            imIA2 = imIA2.*W;
            % compute rms for normalization of cross-correlation
            auxStd1 = stdfast(imIA1);
            auxStd2 = stdfast(imIA2);
            % do the cross-correlation and normalize it
            switch lower(pivPar.ccMethod)
                case 'fft'
                    cc = fftshift(real(ifft2(conj(fft2(imIA1)).*fft2(imIA2))))/(auxStd1*auxStd2)/(iaSizeX*iaSizeY);
                    % find the cross-correlation peak
                    [auxPeak,Upx] = max(max(cc));
                    [aux,Vpx] = max(cc(:,Upx));     %#ok<ASGLU>
                case 'dcn'
                    cc = dcn(imIA1,imIA2,pivPar.ccMaxDCNdist)/(auxStd1*auxStd2)/(iaSizeX*iaSizeY);
                    % find the cross-correlation peak
                    [auxPeak,Upx] = max(max(cc));
                    [aux,Vpx] = max(cc(:,Upx));     %#ok<ASGLU>
                    if (Upx~=iaSizeX/2+ccPxShiftX) || (Vpx~=iaSizeY/2+ccPxShiftY)
                        cc = fftshift(real(ifft2(conj(fft2(imIA1)).*fft2(imIA2))))/(auxStd1*auxStd2)/(iaSizeX*iaSizeY);
                        % find the cross-correlation peak
                        [auxPeak,Upx] = max(max(cc));
                        [aux,Vpx] = max(cc(:,Upx));     %#ok<ASGLU>
                    end
            end
            
            % if the displacement is too large (too close to border), set fail flag
            if (abs(Upx-iaSizeX/2-ccPxShiftX) > pivPar.ccMaxDisplacement*iaSizeX) || ...
                    (abs(Vpx-iaSizeY/2-ccPxShiftY) > pivPar.ccMaxDisplacement*iaSizeY)
                failFlag =  bitset(failFlag,2);
            end
            % corect cc peak for bias caused by interrogation window (see ref. [1], p. 356, eq. (8.104))
            if pivPar.ccCorrectWindowBias && ~isnan(F)
                ccCor = cc./F;
            else
                ccCor = cc;
            end
               % note: this correction is applied only before finding peak position, otherwise spurious peaks are found at
               % borders of IA
            % sub-pixel interpolation (2x3point Gaussian fit, eq. 8.163, p. 375 in [1])
            try
                dU = (log(ccCor(Vpx,Upx-1)) - log(ccCor(Vpx,Upx+1)))/...
                    (log(ccCor(Vpx,Upx-1))+log(ccCor(Vpx,Upx+1))-2*log(ccCor(Vpx,Upx)))/2;
                dV = (log(ccCor(Vpx-1,Upx)) - log(ccCor(Vpx+1,Upx)))/...
                    (log(ccCor(Vpx-1,Upx))+log(ccCor(Vpx+1,Upx))-2*log(ccCor(Vpx,Upx)))/2;
            catch     %#ok<*CTCH>
                failFlag = bitset(failFlag,3);
                dU = NaN; dV = NaN;
            end
            % if imaginary, set fail flag
            if (~isreal(dU)) || (~isreal(dV))
                failFlag = bitset(failFlag,3);
            end
        else
            cc = zeros(iaSizeY,iaSizeX) + NaN;
            auxPeak = NaN;            
            auxStd1 = NaN;            
            auxStd2 = NaN;
            auxMean1 = NaN;            
            auxMean2 = NaN;
            Upx = iaSizeX/2;
            Vpx = iaSizeY/2;
        end
        % save the pivData information about cross-correlation, rough peak position and peak level
        if failFlag == 0
            U(ky,kx) = pivData.iaU0(ky,kx) + Upx + dU - iaSizeX/2 - ccPxShiftX;               % this is subroutine's output
            V(ky,kx) = pivData.iaV0(ky,kx) + Vpx + dV - iaSizeY/2 - ccPxShiftY;               % this is subroutine's output
        else
            U(ky,kx) = NaN;
            V(ky,kx) = NaN;
        end
        status(ky,kx) = failFlag;
        ccPeakIm(1+(ky-1)*iaSizeY:ky*iaSizeY,1+(kx-1)*iaSizeX:kx*iaSizeX) = cc;
        ccPeak(ky,kx) = auxPeak;
        ccStd1(ky,kx) = auxStd1;
        ccStd2(ky,kx) = auxStd2;
        ccMean1(ky,kx) = auxMean1;
        ccMean2(ky,kx) = auxMean2;
        % find secondary peak
        try
            cc(Vpx-2:Vpx+2,Upx-2:Upx+2) = 0;
            ccPeakSecondary(ky,kx) = max(max(cc));
        catch
            try    
                cc(Vpx-1:Vpx+1,Upx-1:Upx+1) = 0;
                ccPeakSecondary(ky,kx) = max(max(cc));
            catch
                ccPeakSecondary(ky,kx) = NaN;
            end
        end % end of secondary peak search
    end % end of loop for ky
end % end of loop for kx

% get IAs where CC failed, and coordinates of corresponding IAs
ccFailedI = logical(bitget(status,2));
ccSubpxFailedI = logical(bitget(status,3));


%% 3. Output results
pivData.Status = uint16(status);
pivData.U = U;
pivData.V = V;
pivData.ccPeak = ccPeak;
pivData.ccPeakSecondary = ccPeakSecondary;
pivData.ccStd1 = ccStd1;
pivData.ccStd2 = ccStd2;
pivData.ccMean1 = ccMean1;
pivData.ccMean2 = ccMean2;
pivData.ccFailedN = sum(sum(ccFailedI));
pivData.ccSubpxFailedN = sum(sum(ccSubpxFailedI));
pivData.ccW = W;
pivData = rmfield(pivData,'iaU0');
pivData = rmfield(pivData,'iaV0');

end


%% LOCAL FUNCTIONS

function [out] = stdfast(in)
% computes root-mean-square (reprogramed, because std in Matlab is somewhat slow due to some additional tests)
in = reshape(in,1,numel(in));
notnan = ~isnan(in);
n = sum(notnan);
in(~notnan) = 0;
avg = sum(in)/n;
out = sqrt(sum(((in - avg).*notnan).^2)/(n-0)); % there should be -1 in the denominator for true std
end


function [cc] = dcn(X1,X2,MaxD)
% computes cross-correlation using discrete convolution
Nx = size(X1,2);
Ny = size(X1,1);
cc = zeros(Ny,Nx);
% create variables defining where is cc(0,0)
dx0 = Nx/2;
dy0 = Ny/2;
if rem(Nx,2) == 0
    dx0 = dx0+1;
else
    dx0 = dx0+0.5;
end
if rem(Ny,2) == 0
    dy0 = dy0+1;
else
    dy0 = dy0+0.5;
end
% pad IAs
X1p = zeros(Ny+2*MaxD,Nx+2*MaxD);
X2p = zeros(Ny+2*MaxD,Nx+2*MaxD);
X1p(MaxD+1:MaxD+Ny,MaxD+1:MaxD+Nx) = X1;
X2p(MaxD+1:MaxD+Ny,MaxD+1:MaxD+Nx) = X2;
% convolve
for kx = -MaxD:MaxD
    for ky = -MaxD:MaxD
        if abs(kx)+abs(ky)>MaxD, continue; end
        cc(dy0+ky,dx0+kx) = sum(sum(...
            X2p(ky+MaxD+1 : ky+MaxD+Ny,  kx+MaxD+1 : kx+MaxD+Nx) .* ...
            X1p(   MaxD+1 : MaxD+Ny,        MaxD+1 : MaxD+Nx)));
    end
end

end

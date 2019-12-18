function [exIm1,exIm2,pivData] = pivInterrogate(im1,im2,pivData,pivPar)
% pivInterrogate - splits images into interrogation areas and create expanded images suitable for pivCrossCorr
%
% Usage:
% [exIm1,exIm2,pivData] = pivInterrogate(im1,im2,pivData,pivPar,Xest,Yest,Uest,Vest)
%
% Inputs:
%    im1,im2 ... image pair (either images, or strings containing paths to image files)
%    pivData ... (struct) structure containing detailed results. No fields are required as the input, but
%             the existing fields will be copies to pivData at the output. If exist, folowing fields will be
%             used:
%        X, Y ... position, at which U, V velocity/displacement is provided 
%        U, V ... displacements in x and y direction (will be used for image deformation). If these fields do
%            not exist, zero velocity is assumed
%    pivPar ... (struct) parameters defining the evaluation. Following fields are considered:
%       iaSizeX, iaSizeY ... size of interrogation area [px]
%       iaStepX, iaStepY ... step between interrogation areas [px]
%       imMask1, imMask2 ... Masking images for im1 and im2. It should be either empty (no mask), or of the
%                            same size as im1 and im2. Masked pixels should be 0 in .imMaskX, non-masked 
%                            pixels should be 1.
%       iaMethod ... way, how interrogation area are created. Possible values are
%           'basic' ... interrogatio areas are regularly distribute rectangles 
%           'offset' ... (not coded) interrogation areas are shifted by the estimated displacement
%           'deflinear' ... (not coded) deformable interrogation areas with linear deformation
%           'defspline' ... (not coded) deformable interrogation areas with spline deformation
%           Note: if Uest and Vest contains only zeros or if they are empty/unspecified, 'basic' method is 
%                 always invoked regardless .iaMethod setting
%       iaImageToDeform ... defines, which image should deform (if .iaMethod is 'deflinear' or 'defspline'), or
%                           in which images IA should be shifted (.iaMethod == 'offset'). Possible values are
%           'image1', 'image2' ... either im1 or im2 is deformed correspondingly to Uest and Vest
%           'both' ... deformation both images are deformed by Uest/2 and Vest/2. More CPU time is required.
%       iaImageInterpolationMethod ... way, how the images are interpolated when deformable IAs are used (for 
%                                      .iaMethod == 'deflinear' or 'defspline'. Possible values are:
%           'linear', 'spline' ... interpolation is carried out using interp2 function with option either
%                                  '*linear' or '*spline'
%       iaPreprocMethod ... defines image preprocessing method. Possible values are
%            'none' ... no image preprocessing
%            'MinMax' ... MinMax filter is applied (see p. 248 in Ref. [1])
%       iaMinMaxSize ... (applies only if iaPreprocMethod is 'MinMax'). Size of MinMax filter kernel.
%       iaMinMaxLevel ... (applies only if iaPreprocMethod is 'MinMax'). Contrast level, below which
%            contrast in not more enhanced.
% Outputs:
%    exIm1, exIm2 ... expanded image 1 and 2. Expanded image is an image, in which IAs are side-by-side 
%          (expanded image has size [iaNX*iaSizeX,iaNY*iaSizeY]); if a pixel appears in n IAs, it will be 
%          present n times in the expanded image. If iaStepX == iaSizeX and iaStepY == iaSizeY and method is 
%          'basic' (or Uest == Vest ==0), expanded image is the same as original image (except cropping).
%    pivData  ... (struct) structure containing more detailed results. If some fiels were present in pivData at the
%              input, they are repeated. Followinf fields are added:
%        imFilename1, imFilename2 ... path and filename of image files (stored only if im1 and im2 are 
%              filenames)
%        imMaskFilename1, imMaskFilename2 ... path and filename of masking files (stored only if imMask1 and 
%              imMask2 are filenames)
%        imArray1, imArray2 ... arrays with read and preprocessed images (these fields are removed at
%              pivAnalyzeImagePair.m)
%        imMaskArray1, imMaskArray2 ... arrays containing Boolean variable (they are removed by
%              pivAnalyzeImagePair.m)
%        imNo1, imNo2, imPairNo ... image number (completed only if im1 and im2 are filenames with images).
%              For example, if im1 and im2 are 'Img000005.bmp' and 'Img000006.bmp', value will be imNo1 = 5, 
%              imNo2 = 6, and imPairNo = 5.5.
%        N ... number of interrogation area (= of velocity vectors)
%        X, Y ... matrices with centers of interrogation areas
%        Status ... matrix with statuis of velocity vectors (uint8). Bits have this coding:
%            1 ... masked (set by pivInterrogate)
%            2 ... cross-correlation failed (set by pivCrossCorr)
%            4 ... peak detection failed (set by pivCrossCorr)
%            8 ... indicated as spurious by median test (set by pivValidate)
%           16 ... interpolated (set by pivReplaced)
%           32 ... smoothed (set by pivSmooth)
%        iaU0, iaV0 ... mean shift of IAs in the deformed image
%        iaSizeX, iaSizeY, iaStepX, iaStepY ... copy of dorresponding fields in pivPar input
%        imSizeX, imSizeY ... image size in pixels
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
%
%% Acronyms and meaning of variables used in this subroutine:
%    IA ... concerns "Interrogation Area"
%    im ... image
%    dx ... some index
%    ex ... expanded (image)
%    est ... estimate (velocity from previous pass) - will be used to deform image
%    aux ... auxiliary variable (which is of no use just a few lines below)
%    cc ... cross-correlation
%    vl ... validation
%    sm ... smoothing
%    Word "velocity" should be understood as "displacement"



%% 0. Read images, if required. Read mask images. Preprocess images, if required
% Extract from pivPar some frequently used fields (for shortenning the code);
iaSizeX = pivPar.iaSizeX;
iaSizeY = pivPar.iaSizeY;
iaStepX = pivPar.iaStepX;
iaStepY = pivPar.iaStepY;

% read images if im1 and im2 are filepaths
if ~isfield(pivData,'imArray1')   % read files only if not read in the previous pass
    if ischar(im1)
        [imgNo] = treatImgPath(im1);
        pivData.imFilename1 = im1;
        pivData.imNo1 = imgNo;
        im1 = imread(im1);
    end
    if ischar(im2)
        [imgNo] = treatImgPath(im2);
        pivData.imFilename2 = im2;
        pivData.imNo2 = imgNo;
        im2 = imread(im2);
    end
    try
        pivData.imPairNo = (pivData.imNo1 + pivData.imNo2)/2;
    catch      %#ok<CTCH>
    end
    
    im1 = single(im1);
    im2 = single(im2);
    
    % read image masks if pivPar.imMask1 and pivPar.imMask2 are filepaths
    if ischar(pivPar.imMask1)
        pivData.imMaskFilename1 = pivPar.imMask1;
        if ~isempty(pivPar.imMask1), imMask1 = imread(pivPar.imMask1);
        else
            imMask1 = [];
        end
    else
        imMask1 = pivPar.imMask1;
    end
    if ischar(pivPar.imMask2)
        pivData.imMaskFilename2 = pivPar.imMask2;
        if ~isempty(pivPar.imMask2), imMask2 = imread(pivPar.imMask2);
        else
            imMask2 = [];
        end
    else
        imMask2 = pivPar.imMask2;
    end
    pivData.imMaskArray1 = imMask1;
    pivData.imMaskArray2 = imMask2;
    
    % check the consistence of images
    auxDiff = abs(size(im1)-size(im2))+abs(size(imMask1)-size(imMask2));
    if numel(imMask1) > 0
        auxDiff = auxDiff + abs(size(im1)-size(imMask1));
    end
    if sum(auxDiff)>0
        error('PIVsuite:InconsistImgs','Image 1, Image 2 (and possible mask image) are inconsistent in size.');
    end
    
    % Preprocess images.
    if strcmpi(pivPar.iaPreprocMethod,'minmax')
        im1 = minMaxFilter(im1,pivPar,imMask1);
        im2 = minMaxFilter(im2,pivPar,imMask2);
    end
    pivData.imArray1 = im1;
    pivData.imArray2 = im2;
    
else
    % if images read and preprocessed in previous pass, read them from pivData variable
    im1 = pivData.imArray1;
    im2 = pivData.imArray2;
    imMask1 = pivData.imMaskArray1;
    imMask2 = pivData.imMaskArray2;
end

%% 1. Compute position of IA's
% get size of the image
imSizeX = size(im1,2);
imSizeY = size(im1,1);

% get the number of IA's
iaNX = floor((imSizeX - iaSizeX)/iaStepX)+1;
iaNY = floor((imSizeY - iaSizeY)/iaStepY)+1;

% distribute IA's (undeformed image, no offset):
auxLengthX = iaStepX * (iaNX-1) + iaSizeX;
auxLengthY = iaStepY * (iaNY-1) + iaSizeY;
auxFirstIAX = floor((imSizeX - auxLengthX)/2) + 1;
auxFirstIAY = floor((imSizeY - auxLengthY)/2) + 1;
iaStartX = (auxFirstIAX:iaStepX:(auxFirstIAX+(iaNX-1)*iaStepX))';  % first columns of IA's
iaStartY = (auxFirstIAY:iaStepY:(auxFirstIAY+(iaNY-1)*iaStepY))';  % first rows of IA's
iaStopX = iaStartX + iaSizeX - 1;  % last columns of IA's
iaStopY = iaStartY + iaSizeY - 1;  % last rows of IA's
iaCenterX = (iaStartX + iaStopX)/2;    % center of IA's (usually between pixels)
iaCenterY = (iaStartY + iaStopY)/2;
[X,Y] = meshgrid(iaCenterX,iaCenterY); % this is a mesh, at which velocity is detrmined
        % in last line, there is iaCenterX-1, because pixel im(1,1) has coordinates (0,0)

% initialize status variable
status = zeros(iaNY,iaNX);

% if velocity estimate is not specified, initialize it. If velocity field is specified,use it as velocity
% estimation:
if ~isfield(pivData,'X')||~isfield(pivData,'Y')||~isfield(pivData,'U')||~isfield(pivData,'V')
    Xest = X;
    Yest = Y; 
    Uest = zeros(iaNY,iaNX); 
    Vest = zeros(iaNY,iaNX); 
else
    Xest = pivData.X;
    Yest = pivData.Y;
    Uest = pivData.U;
    Vest = pivData.V;
end
% if velocity estimation is zero (happens also if not spcified), set method as 'basic'
if max(max(abs(Uest)))+max(max(abs(Vest)))<20*eps
    pivPar.iaMethod = 'basic';
end
if strcmpi(pivPar.iaMethod,'basic')
    Uest = 0*Uest;
    Vest = 0*Vest;
end


%% 2. Mask images
% mark masked pixels as NaNs. Later, NaNs are replaced by mean of non-masked pixels within each IA. 
if numel(imMask1) > 0
    imMask1 = ~logical(imMask1); 
    im1(imMask1) = NaN;
end
if numel(imMask2) > 0
    imMask2 = ~logical(imMask2); 
    im2(imMask2) = NaN;
end


%% 3. Create expanded images
% Expanded image is an image, in which IAs are side-by-side (expanded image has size 
% [iaNX*iaSizeX,iaNY*iaSizeY]); if a pixel appears in n IAs, it will be present n times in the expanded image.
% If iaStepX == iaSizeX and iaStepY == iaSizeY and method is 'basic' (or Uest == Vest ==0), expanded image is
% the same as original image (except cropping)

% initialize expanded images
exIm1 = single(zeros(iaNY*iaSizeY, iaNX*iaSizeX) + NaN);
exIm2 = exIm1;

% do everything with first image, then with the second

switch lower(pivPar.iaMethod)
    case 'basic'
        % standard interrogation - no IA offset or deformation
        for kx = 1:iaNX
            for ky = 1:iaNY
                % get the interrogation areas
                imIA1 = im1(iaStartY(ky):iaStopY(ky),iaStartX(kx):iaStopX(kx));
                imIA2 = im2(iaStartY(ky):iaStopY(ky),iaStartX(kx):iaStopX(kx));
                % set the masked pixel to mean value of remaining pixels in the IA
                masked1 = isnan(imIA1);
                masked2 = isnan(imIA2);
                auxMean1 = sum(sum(imIA1.*(~masked1)))/sum(sum(~masked1)); % this is man of pixels inside - faster than "mean" function
                auxMean2 = sum(sum(imIA2.*(~masked2)))/sum(sum(~masked2));
                imIA1(masked1) = auxMean1;
                imIA2(masked2) = auxMean2;
                % copy it to expanded images
                exIm1(1+(ky-1)*iaSizeY:ky*iaSizeY,1+(kx-1)*iaSizeX:kx*iaSizeX) = imIA1;
                exIm2(1+(ky-1)*iaSizeY:ky*iaSizeY,1+(kx-1)*iaSizeX:kx*iaSizeX) = imIA2;
                % check the number of masked pixels, and if larget than 1/2*iaSizeX*iaSizeY, % consider IA as masked or outside
                if sum(sum(logical(masked1+masked2))) > 0.5*iaSizeX*iaSizeY;
                    status(ky,kx) = 1;
                end
            end
        end
        % set the interpolated velocity to zeros
        U0 = zeros(iaNY,iaNX);
        V0 = zeros(iaNY,iaNX);
    case 'offset'
        % interrogation with offset of IA, no deformation of IA        
        % interpolate the velocity estimates to the new grid and round it
        if sum(sum(isnan(Uest+Vest)))>0
            Uest = inpaint_nans(Uest);
            Vest = inpaint_nans(Vest);
        end
        U0 = interp2(Xest,Yest,Uest,X,Y,'linear');
        V0 = interp2(Xest,Yest,Vest,X,Y,'linear');
        if sum(sum(isnan(U0+V0)))>0
            U0 = inpaint_nans(U0);
            V0 = inpaint_nans(V0);
        end
        if ~strcmpi(pivPar.iaImageToDeform,'both')
            U0 = round(U0);
            V0 = round(V0);
        else
            U0 = 2*round(U0/2);
            V0 = 2*round(V0/2);
        end
        % create index matrices:  in which position of corresponding pixels in image pair is stored
        [auxX,auxY] = meshgrid(0:iaSizeX-1,0:iaSizeY-1);
        for kx = 1:iaNX
            for ky = 1:iaNY
                % calculate the shift of IAs
                switch lower(pivPar.iaImageToDeform)
                    case 'image1'
                        dxX1 = iaStartX(kx) + auxX - U0(ky,kx);
                        dxX2 = iaStartX(kx) + auxX;
                        dxY1 = iaStartY(ky) + auxY - V0(ky,kx);
                        dxY2 = iaStartY(ky) + auxY;
                    case 'image2'
                        dxX1 = iaStartX(kx) + auxX;
                        dxX2 = iaStartX(kx) + auxX + U0(ky,kx)/2;
                        dxY1 = iaStartY(ky) + auxY;
                        dxY2 = iaStartY(ky) + auxY + V0(ky,kx)/2;
                    case 'both'
                        dxX1 = iaStartX(kx) + auxX - U0(ky,kx)/2;
                        dxX2 = iaStartX(kx) + auxX + U0(ky,kx)/2;
                        dxY1 = iaStartY(ky) + auxY - V0(ky,kx)/2;
                        dxY2 = iaStartY(ky) + auxY + V0(ky,kx)/2;
                end
                % check, where the shifted pixel is out of image (shifted IA goes outside image), and set the
                % corresponding index to 1 (will be corrected later). Pixels outside will be treated as masked
                masked1 = logical(logical(dxX1<1) + logical(dxX1>size(im1,2)) + ...
                    logical(dxY1<1) + logical(dxY1>size(im1,1)));
                masked2 = logical(logical(dxX2<1) + logical(dxX2>size(im2,2)) + ...
                    logical(dxY2<1) + logical(dxY2>size(im2,1)));
                dxX1(masked1) = 1;
                dxY1(masked1) = 1;
                dxX2(masked2) = 1;
                dxY2(masked2) = 1;
                % convert double-indexing of matrix (M(k,j)) to single index (M(l), l = (j-1)*rows + k
                dxS1 = (dxX1-1)*size(im1,1) + dxY1;
                dxS2 = (dxX2-1)*size(im2,1) + dxY2;
                % copy the IA from the original image
                imIA1 = im1(dxS1);
                imIA2 = im2(dxS2);
                % add masked pixel (NaNs) to outside pixels (masked1 and masked2)
                masked1 = logical(masked1 + isnan(imIA1));
                masked2 = logical(masked2 + isnan(imIA2));
                % correct masked pixels - replace them by mean
                auxMean1 = sum(sum(imIA1.*(~masked1)))/sum(sum(~masked1)); % this is mean of pixels inside - faster than "mean" function
                auxMean2 = sum(sum(imIA2.*(~masked2)))/sum(sum(~masked2));
                imIA1(masked1) = auxMean1;
                imIA2(masked2) = auxMean2;
                % copy IA to the expanded image
                exIm1(1+(ky-1)*iaSizeY:ky*iaSizeY,1+(kx-1)*iaSizeX:kx*iaSizeX) = imIA1;
                exIm2(1+(ky-1)*iaSizeY:ky*iaSizeY,1+(kx-1)*iaSizeX:kx*iaSizeX) = imIA2;
                % check the number of masked or outside pixels, and if larget than 1/2*iaSizeX*iaSizeY,
                % consider IA as masked or outside
                if sum(sum(logical(masked1+masked2))) > 0.5*iaSizeX*iaSizeY;
                    status(ky,kx) = 1;
                end
            end
        end
    case {'deflinear','defspline'}
        % cases with deformable IAs
        % create new X and Y mesh. Make it larger (used for extrapolation of velocity estimates for methods
        % 'deflinear' and 'defspline')
        auxCenterX = [iaCenterX(1)-(ceil(iaCenterX(1)/iaStepX):-1:1)'*iaStepX ;...
            iaCenterX ; ...
            iaCenterX(end)+ (1:ceil((imSizeX-iaCenterX(end))/iaStepX))'*iaStepX];
        auxCenterY = [iaCenterY(1)-(ceil(iaCenterY(1)/iaStepY):-1:1)'*iaStepY ;...
            iaCenterY ; ...
            iaCenterY(end)+ (1:ceil((imSizeY-iaCenterY(end))/iaStepY))'*iaStepY];
        [Xextrap,Yextrap] = meshgrid(auxCenterX,auxCenterY);
        % extrapolate velocity estimates to the extrapolation mesh
        auxOK = logical(~isnan(Uest+Vest));
        Uest(~auxOK) = NaN;
        Vest(~auxOK) = NaN;
        Uest = inpaint_nans(Uest);
        Vest = inpaint_nans(Vest);
        try
            % this will run only in new versions of Matlab
            Uestimator = griddedInterpolant(Xest',Yest',Uest','spline');
            Vestimator = griddedInterpolant(Xest',Yest',Vest','spline');
            Uextrap = Uestimator(Xextrap',Yextrap');
            Vextrap = Vestimator(Xextrap',Yextrap');
            Uextrap = inpaint_nans(Uextrap');
            Vextrap = inpaint_nans(Vextrap');
        catch err
            % if error (griddedInterpolant missing in older versions), use TriScatteredInterp
            if (strcmp(err.identifier,'MATLAB:UndefinedFunction'))
                auxXest = reshape(Xest,numel(Xest),1);
                auxYest = reshape(Yest,numel(Yest),1);
                auxUest = reshape(Uest,numel(Uest),1);
                auxVest = reshape(Vest,numel(Vest),1);
                Uestimator = TriScatteredInterp(auxXest,auxYest,auxUest,'natural');   %#ok<DTRIINT>
                Vestimator = TriScatteredInterp(auxXest,auxYest,auxVest,'natural'); %#ok<DTRIINT>
                Uextrap = Uestimator(Xextrap,Yextrap);
                Vextrap = Vestimator(Xextrap,Yextrap);
                Uextrap = inpaint_nans(Uextrap);
                Vextrap = inpaint_nans(Vextrap);
            else
                retrow(err);
            end
        end
        % initialize matrices with subtracted deformation
        U0 = X + NaN;
        V0 = U0;
        % create marices with pixel coordines of deformed image
        % start with undeformed coordinates
        [coordX, coordY] = meshgrid(1:size(im1,2),1:size(im1,1));
        % estimate velocities on these coordinates (one can use also TriScatteredInterp, but interp2 is faster)
        switch lower(pivPar.iaMethod)
            case 'deflinear'
                Udef = interp2(Xextrap,Yextrap,Uextrap,coordX,coordY,'*linear');
                Vdef = interp2(Xextrap,Yextrap,Vextrap,coordX,coordY,'*linear');
            case 'defspline'
                Udef = interp2(Xextrap,Yextrap,Uextrap,coordX,coordY,'*spline');
                Vdef = interp2(Xextrap,Yextrap,Vextrap,coordX,coordY,'*spline');
            otherwise
                error('PIVsuite:iaMethod','Unknown iaMethod.');
        end
        % deform the coordinates
        switch lower(pivPar.iaImageToDeform)
            case 'image1'
                coordX1 = coordX - Udef;
                coordY1 = coordY - Vdef;
                coordX2 = coordX;
                coordY2 = coordY;
                auxIm1Def = true; auxIm2Def = false;   % flags if the corresponding image deforms
            case 'image2'
                coordX1 = coordX;
                coordY1 = coordY;
                coordX2 = coordX + Udef;
                coordY2 = coordY + Vdef;
                auxIm1Def = false; auxIm2Def = true;
            case 'both'
                coordX1 = coordX - Udef/2;
                coordY1 = coordY - Vdef/2;
                coordX2 = coordX + Udef/2;
                coordY2 = coordY + Vdef/2;
                auxIm1Def = true; auxIm2Def = true;
        end
        % find pixels, which are outside real images - will be treated as masked pixels
        masked1 = logical(logical(coordX1<1) + logical(coordX1>size(im1,2)) + ...
            logical(coordY1<1) + logical(coordY1>size(im1,1)));
        masked2 = logical(logical(coordX2<1) + logical(coordX2>size(im2,2)) + ...
            logical(coordY2<1) + logical(coordY2>size(im2,1)));
        % bring pixel outside to inside...
        coordX1(masked1) = 1;
        coordY1(masked1) = 1;
        coordX2(masked2) = 1;
        coordY2(masked2) = 1;
        % deform the image
        switch lower(pivPar.iaImageInterpolationMethod)
            case 'linear'
                if auxIm1Def, im1 = exp(interp2(coordX,coordY,log(im1+1),coordX1,coordY1,'*linear'))-1; end
                if auxIm2Def, im2 = exp(interp2(coordX,coordY,log(im2+1),coordX2,coordY2,'*linear'))-1; end
            case 'spline'
                if auxIm1Def   % interpolate images only if they deform - it is CPU consuming
                    % mask is lost - NaNs in images disappear with interp2 is method is '*spline'! :-(
                    im1(isnan(im1)) = 0;
                    im1 = exp(interp2(coordX,coordY,log(im1+1),coordX1,coordY1,'*spline'))-1;
                    % Masking images are therefore interpolated separately:
                    if numel(imMask1)>0
                        auxMaskDef = interp2(coordX,coordY,single(imMask1),coordX1,coordY1,'*spline');
                        auxMaskDef(isnan(auxMaskDef)) = 0;
                        auxMaskDef = logical(auxMaskDef>0.5);
                        im1(auxMaskDef) = NaN;
                    end
                end % end of "image 1 deforms"
                if auxIm2Def
                    im2(isnan(im2)) = 0;
                    im2 = exp(interp2(coordX,coordY,log(im2+1),coordX2,coordY2,'*spline'))-1;
                    if numel(pivPar.imMask2)>0
                        auxMaskDef = interp2(coordX,coordY,single(imMask2),coordX2,coordY2,'*spline');
                        auxMaskDef(isnan(auxMaskDef)) = 0;
                        auxMaskDef = logical(auxMaskDef>0.5);
                        im1(auxMaskDef) = NaN;
                    end
                end % end of "Image 2 deforms"
            otherwise
                error('PIVsuite:iaImageInterpolationMethod','Unkonwn iaImageInterpolationMethod');
        end
        % set pixels outside image to NaN
        im1(masked1) = NaN;
        im2(masked2) = NaN;
        % create the expanded images
        for kx = 1:iaNX
            for ky = 1:iaNY
                % take IA
                imIA1 = im1(iaStartY(ky):iaStopY(ky),iaStartX(kx):iaStopX(kx));
                imIA2 = im2(iaStartY(ky):iaStopY(ky),iaStartX(kx):iaStopX(kx));
                % replace pixels with NaNs (either masked or outside image) by mean of remaining pixels
                masked1 = isnan(imIA1);
                masked2 = isnan(imIA2);
                imIA1(masked1) = 0;
                imIA2(masked2) = 0;
                auxMean1 = sum(sum(imIA1))/sum(sum(~masked1)); % this is mean of non-NaN pixels - it is faster than mean.m
                auxMean2 = sum(sum(imIA2))/sum(sum(~masked2));
                imIA1(masked1) = auxMean1;
                imIA2(masked2) = auxMean2;
                % store in the expanded image
                exIm1(1+(ky-1)*iaSizeY:ky*iaSizeY,1+(kx-1)*iaSizeX:kx*iaSizeX) = imIA1;
                exIm2(1+(ky-1)*iaSizeY:ky*iaSizeY,1+(kx-1)*iaSizeX:kx*iaSizeX) = imIA2;
                % compute the mean shift of deformed IAs
                auxU = Udef(iaStartY(ky):iaStopY(ky),iaStartX(kx):iaStopX(kx));
                auxV = Vdef(iaStartY(ky):iaStopY(ky),iaStartX(kx):iaStopX(kx));
                U0(ky,kx) = sum(sum(auxU.*logical((~masked1).*(~masked2))))/sum(sum(logical((~masked1).*(~masked2))));
                V0(ky,kx) = sum(sum(auxV.*logical((~masked1).*(~masked2))))/sum(sum(logical((~masked1).*(~masked2))));
                % consider IA as masked if its center is masked in both images
                auxM1 = indexi(imMask1,Y(ky,kx),X(ky,kx));
                auxM2 = indexi(imMask2,Y(ky,kx),X(ky,kx));
                if auxM1+auxM2>=1 ||...
                ...% check the number of masked or border pixels, and if larger than 1/3*iaSizeX*iaSizeY,
                ...% consider IA as masked 
                    sum(sum(logical(masked1+masked2))) > 0.5*iaSizeX*iaSizeY;
                    status(ky,kx) = 1;
                end
            end
        end
end


%% 4. Output results via pivData variable
pivData.X = X;
pivData.Y = Y;
pivData.U = X + NaN;
pivData.V = X + NaN;
pivData.N = numel(X);
pivData.Status = uint16(status);
pivData.maskedN = sum(sum(logical(bitget(status,1))));
pivData.imSizeX = imSizeX;
pivData.imSizeY = imSizeY;
pivData.iaSizeX = iaSizeX;
pivData.iaSizeY = iaSizeY;
pivData.iaStepX = iaStepX;
pivData.iaStepY = iaStepY;
pivData.iaU0 = U0;
pivData.iaV0 = V0;
end


%% LOCAL FUNCTIONS

function [ corrected ] = minMaxFilter(im,pivPar,Mask)
% MinMax filter - corrects uneven background and normalizes image contrast
%   adapted following algorithm described on p. 250, Ref. [1]
%
%   This subroutine requires image processing toolbox
%
S = pivPar.iaMinMaxSize;
L = pivPar.iaMinMaxLevel;
% create masking matrix (ones in a circular matrix)
domain = ones(S,S);
auxX = ones(S,1)*(-(S-1)/2:1:(S-1)/2);
auxY = (-(S-1)/2:1:(S-1)/2)'*ones(1,S);
auxD = sqrt(auxX.^2+auxY.^2);
domain(auxD+1/4 >= (S-1)/2) = 0;
N = sum(sum(domain));
domain = single(domain);
% if no Mask is specified, skip it
if nargin<3
    Masking = false;
elseif numel(Mask)==0
    Masking = false;
else
    % check if masking is present
    Mask = (Mask==0);    % Masked pixel are 1 now
    Masking = sum(sum(Mask))>0;
end
% Compute local low value and filter it
im = double(im);
if Masking
    im1 = im;
    im1(Mask) = max(max(im1));
    Lo = ordfilt2(im1,1,domain,'symmetric');
    Lo = filter2(domain,Lo)/N;
else
    Lo = filter2(domain,ordfilt2(im,1,domain,'symmetric'))/N;
end
% Compute local high value and filter it
if Masking
    im1 = im;
    im1(Mask) = min(min(im1));
    Hi = ordfilt2(im1,N,domain,'symmetric');
    Hi = filter2(domain,Hi)/N;
else
    Hi = filter2(domain,ordfilt2(im,N,domain,'symmetric'))/N;
end
% enlarge mask (pixels in enlarged mask will not be considered during normalization)
if Masking
    MaskF = imdilate(Mask,domain);
end
% compute contrast and put lower limit on it
contrast = Hi-Lo;
contrast = (contrast>L).*(contrast-L) + L;
corrected = single((im-Lo)./contrast);
% normalize image
corrMax = corrected;
if Masking, corrMax(MaskF) = 0;end
corrMax = max(max(corrMax(S+1:end-S,S+1:end-S)));
corrected = 255*corrected/corrMax;
corrected(corrected>255) = single(255);
end


function [imgNo, filename, folder] = treatImgPath(path)
% separate the path to get the folder, filename, and number if contained in the name
filename = '';
imgNo = [];
folder = '';
if numel(path)>0
    path = path(end:-1:1);
    I = find(path=='/'|path=='\');
    I = I(1);
    Idot = find(path=='.');
    Idot = Idot(1);
    try
        folder = path(I+1:end);
        folder = folder(end:-1:1);
    catch  %#ok<CTCH>
        folder = '';
    end
    try
        filename = path(Idot+1:I-1);
        filename = filename(end:-1:1);
    catch  %#ok<CTCH>
        filename = '';
    end
    try
        aux = regexp(filename,'[0-9]');
        aux = filename(aux);
        imgNo = str2double(aux);
    catch  %#ok<CTCH>
        imgNo = [];
    end
end
end


function [out] = indexi(array,I1,I2)
% index array for getting values between items (if In is not integer, give value on midway betweem pixels);
if numel(array) == 0
    out = 0; 
    return;
end
A1 = array(floor(I1),floor(I2));
A2 = array(ceil(I1),floor(I2));
A3 = array(floor(I1),ceil(I2));
A4 = array(ceil(I1),ceil(I2));
out = (A1+A2+A3+A4)/4;
end

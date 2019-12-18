function [pivData,ccFunction] = pivAnalyzeImagePair(im1,im2,pivPar)
% pivAnalyzeImagePair - performs single- or multi-pass analysis of displacement between two images using PIV technique
%
% Usage:
% [pivData] = pivAnalyzeImagePair(im1,im2,pivData,pivDataIn);
% [pivData,ccFunction] = = pivAnalyzeImagePair(im1,im2,pivData,pivDataIn);
%
% Inputs:
%    im1,im2 ... image pair (either images, or strings containing paths to image files)
%    pivDataIn ... (struct) structure containing detailed results. No fields are required as the input, but
%          the existing fields will be copies to pivData at the output. If exist, folowing fields will be used:
%        X, Y ... position, at which U, V velocity/displacement is provided
%        U, V ... displacements in x and y direction (will be used for image deformation). If these fields do
%            not exist, zero velocity is assumed
%    pivPar ... (struct) parameters defining the evaluation. Use pivParams.m for creating pivPar. Following
%            fields are considered:
%      --- prefix an - these fields are used mostly in subroutines pivAnalyseImagePair.m
%          anNpasses ... number of passes
%      --- prefix ia - these fields are used mostly in subroutine pivInterrogate.m
%          iaSizeX, iaSizeY ... size of interrogation area [px]
%          iaStepX, iaStepY ... step between interrogation areas [px]
%          imMask1, imMask2 ... Masking images for im1 and im2. It should be either empty (no mask), or of the
%               same size as im1 and im2. Masked pixels should be 0 in .imMaskX, non-masked pixels should be 1
%          iaMethod ... way, how interrogation area are created. Possible values are
%               'basic' ... interrogatio areas are regularly distribute rectangles
%               'offset' ... (not coded) interrogation areas are shifted by the estimated displacement
%               'deflinear' ... (not coded) deformable interrogation areas with linear deformation
%               'defspline' ... (not coded) deformable interrogation areas with spline deformation
%             - Note: if Uest and Vest contains only zeros or if they are empty/unspecified, 'basic' method is
%                    always invoked regardless .iaMethod setting
%          iaImageToDeform ... defines, which image should deform. It is taken in account if .iaMethod is
%                  'deflinear' or 'defspline', or if it is 'offset' (then it defines, in which image IAs are
%                  shifted). Possible values are
%              'image1', 'image2' ... either im1 or im2 is deformed correspondingly to Uest and Vest
%              'both' ... deformation both images are deformed by Uest/2 and Vest/2. Sligthly more CPU time is
%                         required.
%          iaImageInterpolationMethod ... way, how the images are interpolated when deformable IAs are used
%                  (for .iaMethod == 'deflinear' or 'defspline'. Possible values are:
%              'linear', 'spline' ... interpolation is carried out using interp2 function with option either
%                                     '*linear' or '*spline'
%      --- prefix cc - these fields are used mostly in subroutines pivCrossCorr.m
%          ccRemoveIAMean ... if =0, do not remove IA's mean before cross-correlation; if =1, remove the mean;
%              if in between, remove the mean partially
%          ccMaxDisplacement ... maximum allowed displacement to accept cross-correlation peak. This parameter
%              is a multiplier; e.g. setting ccMaxDisplacement = 0.6 means that the cross-correlation peak
%              must be located within [-0.6*iaSizeX...0.6*iaSizeX, -0.6*iaSizeY...0.6*iaSizeY] from the zero
%              displacement.
%            - Note: IA offset is not included in the displacement accounted by ccMaxDisplacement, hence real
%              displacement can be larger if ccIAmethod is any other than 'basic'.
%          ccWindow ... windowing function p. 389 in ref [1]). Possible values are 'uniform', 'gauss',
%             'parzen', 'hanning', 'welch'
%            - default value: 'Welch'. If missing in pivPar during evaluation, 'unifor' is considered for
%              compatibility reasons.
%          ccCorrectWindowBias ... Set whether cc is corrected for bias due to IA windowing.
%            - default value: true (will correct peak)
%      --- prefix vl - these fields relates to vector validation, subroutine pivValidate.m
%          vlTresh, vlEps ... Define treshold for the median test. To accepted, the difference of actual vector
%              from the median (of vectors in the neighborhood) should be at most vlTresh *(vlEps +
%              (neighborhood vectors) - (their median))
%          vlDist ... to what distance median test is performed (if vlDist = 1, kernel has size 3x3; for
%              vlDist = 2, kernel is 5x5, and so on)
%          vlPasses ... number of passes of the median test
%      --- prefix rp - affects replacement of invalid vectors, subroutine pivReplace.m
%          rpMethod ... specifies how the spurious vectors are replaced. Possible values are
%              'none' ... do not replace spurious vectors
%              'linear' ... replace spurious vectors with the use of TriScatteredInterp in Matlab, specifying
%                   its method as 'linear'. If pivData contains data for image sequence, replacement is done 
%                   in each time slice indepedently on other time slices.
%              'natural' ... replace spurious vectors with the use of TriScatteredInterp in Matlab, specifying
%                   its method as 'natural'. If pivData contains data for image sequence, replacement is done 
%                   in each time slice indepedently on other time slices.
%              'inpaint' ... use D'Errico's subroutine "inpaint_nans". If pivData contains data for image 
%                   sequence, replacement is done in each time slice indepedently on other time slices.
%              'inpaintGarcia' ... use Garcia's subroutine "inpaintn". If pivData contains data for image 
%                   sequence, replacement is done in each time slice indepedently on other time slices.
%              'linearT' ... replace spurious vectors with the use of TriScatteredInterp in Matlab, specifying
%                   its method as 'linear'. If pivData contains data for image sequence, replacement considers
%                   also values in other time slices.
%              'naturalT' ... replace spurious vectors with the use of TriScatteredInterp in Matlab, 
%                   specifying its method as 'natural'. If pivData contains data for image sequence, 
%                   replacement considers also values in other time slices.
%              'inpaintT' ... use D'Errico's subroutine "inpaint_nans". If pivData contains data for image 
%                   sequence, replacement considers also values in other time slices.
%              'inpaintGarciaT' ... use Garcia's subroutine "inpaintn". If pivData contains data for image 
%                   sequence, replacement considers also values in other time slices.
%      --- prefix sm - affects smoothing of vector field, subroutine pivSmooth.m
%          smMethod ... defines smoothing method. Possible values are:
%              'none' ... do not perform smoothing
%              'smoothn' ... uses smoothn.m function by Damian Garcia [5]
%              'gauss' ... uses Gaussian kernel
%          smSigma ... amount of smoothing
%          smSize ... size of filter (applies only to Gaussian filter)
%      --- Following fields (with qv prefix) are used to plot results during computing.
%          qvOptionsPair ... (cell array) string of options that is send to pivQuiver.m to plot intermediate
%              results after each processing pass by pivAnalyseImagePair.m
%            - for options for pivQuiver, see help of pivQuiver
%
% Outputs:
%    pivData  ... (struct) structure containing more detailed results. If some fiels were present in pivData
%           at the input, they are repeated. Followinf fields are added:
%        imFilename1, imFilename2 ... path and filename of image files (stored only if im1 and im2 are
%              filenames)
%        imMaskFilename1, imMaskFilename2 ... path and filename of masking files (stored only if imMask1 and
%              imMask2 are filenames)
%        N ... number of interrogation area (= of velocity vectors)
%        X, Y ... matrices with centers of interrogation areas (positions of velocity vectors)
%        U, V ... components of velocity vectors
%        Status ... matrix with statuis of velocity vectors (uint8). Bits have this coding:
%            1 (bit 1) ... masked (set by pivInterrogate)
%            2 (bit 2) ... cross-correlation failed (set by pivCrossCorr)
%            4 (bit 3) ... peak detection failed (set by pivCrossCorr)
%            8 (bit 4) ... indicated as spurious by median test based on image pair (set by pivValidate)
%           16 (bit 5) ... interpolated (set by pivReplaced)
%           32 (bit 6) ... smoothed (set by pivSmooth)
%           64 (bit 7) ... indicated as spurious by median test based on image sequence (set by pivValidate);
%              this flag cannot be set when working with a single image pair
%          128 (bit 8) ... interpolated within image sequence (set by pivReplaced); this flag cannot be set 
%              when working with a single image pair
%          256 (bit 9) ... smoothed within an image sequence (set by pivSmooth); this flag cannot be set when
%              working with a single image pair
%           (example: if Status for a particulat point is 56 = 32 + 16 + 8, the velocity vector in this point
%            was indicated as spurious, was replaced by interpolating neighborhood values and was then
%            adjusted by smoothing.)
%        iaSizeX, iaSizeY, iaStepX, iaStepY ... copy of dorresponding fields in pivPar input
%        imSizeX, imSizeY ... image size in pixels
%        imFilename1, imFilename2 ... path and filename of image files (stored only if im1 and im2 are
%            filenames)
%        imMaskFilename1, imMaskFilename2 ... path and filename of masking files (stored only if imMask1 and
%            imMask2 are filenames)
%        imNo1, imNo2, imPairNo ... image number and number of image pair (stored only if im1 and im2 are
%            string with filenames of images). For example, if im1 and im2 are 'Img000005.bmp' and
%            'Img000006.bmp', value will be imNo1 = 5, imNo2 = 6, and imPairNo = 5.5.
%        ccPeak ... table with values of cross-correlation peak
%        ccPeakSecondary ... table with values of secondary cross-correlation peak (maximum of
%                            crosscorrelation, if 5x5 neighborhood of primary peak is removed)
%        ccFailedN ... number of vectors for which cross-correlation failed
%            at distance larger than ccMaxDisplacement*(iaSizeX,iaSizeY) )
%        ccSubpxFailedN ... number of vectors for which subpixel interpolation failed
%        spuriousN ... number of spurious vectors (status 1)
%        spuriousX, spuriousY ... positions, at which the velocity is spurious
%        spuriousU, spuriousV ... components of the velocity/displacement vectors, which were indicated as
%                             spurious
%        replacedN ... number of interpolated vectors (status 2)
%        replacedX,replacedY ... positions, at which velocity/displacement vectors were replaced
%        replacedU,replacedV ... components of the velocity/displacement vectors, which were replaced
%        validN ... number of original and vectors
%        validX,validY ... positions, at which velocity/displacement vectors is original and valid
%        validU,validV ... original and valid components of the velocity/displacement vector
%        infCompTime ... 1D array containing computational time of individual passes (in seconds)
%    ccFunction ... returns cross-correlation function (in form of an expanded image); see pivCrossCorr
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
%%
pivPar = [];      % variable for settings
pivData = [];     % variable for storing results

pivPar.imMask1 = [];      % No mask for first image
pivPar.imMask2 = [];      % No mask for the second image

%% Image preprocessing

pivPar.iaPreprocMethod = 'MinMax';  % include image preprocessing by local MinMax filter
pivPar.iaMinMaxSize = 7;            % set kernel size of MinMax filter
pivPar.iaMinMaxLevel = 16;          % set minimal level, below which the contrast is not improved 


%% Set the *windowing function*:

pivPar.ccWindow = 'Gauss2';   % This filter is relatively narrow and will 
                              % help measure velocities in middle of IA's
                              % even for large IA'a     

%% Vector validation settings

pivPar.vlMinCC = 0.3;      % reject all vectors with normalized 
                           % cross-correlation smaller than 0.3 of median
                           
pivPar.vlPasses = [1 1 1 1 1];        % number of passes
pivPar.vlDist = 2 ;               % distance, to which vectors are included in the local
                                  % statistics (vlDist = 2 stays for 5x5 vector neighborhood
pivPar.vlTresh = 3;               % threshold value
pivPar.vlEps = 0.08;              % allowed errors

%% Vector smoothing settings

pivPar.smMethod = 'smoothn';    % choose the method (possible values 'none', 'Gauss', 'smoothn')
pivPar.smSigma = 0.2;         % larger value means smoother results; NaN stays for automatic setting
% pivPar.smSize = 10;            % size of the smoothing filter (applies only to 
                              % gaussian filter)
                              
%% Set remaining parameters and run the analysis

[pivParIn, pivData] = pivParams(pivData,pivPar,'defaults');   
pivData.infCompTime = [];
% loop for all required passes
for kp = 1:5
    timer = tic;
    pivData0 = pivData;     % save velocity before computation - will be used if predictor-corrector is used
    % extract parameters for the corresponding pass
    [pivPar] = pivParams(pivData,pivParIn,'singlePass',kp);
    % find interrogation areas in images, shift or deform them if required
    [exIm1,exIm2,pivData] = pivInterrogate(im1,im2,pivData,pivPar);
    % compute cross-correlation between interrogation areas
    if (kp == pivParIn.anNpasses) && nargout > 1
        [pivData,ccFunction] = pivCrossCorr(exIm1,exIm2,pivData,pivPar);
    else
        pivData = pivCrossCorr(exIm1,exIm2,pivData,pivPar);
    end
    % apply predictor-corrector to the velocity data
    pivData = pivCorrector(pivData,pivData0,pivPar);
    % validate velocity field
    pivData = pivValidate(pivData,pivPar);
    % interpolate invalid velocity vectros
    pivData = pivReplace(pivData,pivPar);
    % smooth the velocity field (if smoothing is not required, set pivDataIn.smMethod = 'none')
    pivData = pivSmooth(pivData,pivPar);
    % save the information about actual pass
    pivData.infPassNo = kp;
    % show the plot if reqquired
    if ~(~usejava('jvm') || ~usejava('desktop') || ~feature('ShowFigureWindows'))
        if isfield(pivPar,'qvPair') && numel(pivPar.qvPair)>0
            pivQuiver(pivData,pivPar.qvPair);
            title(['Pass no. ',num2str(kp,'%d'),...
                ', IA size ', num2str(pivPar.iaSizeX,'%d'), ...
                'x', num2str(pivPar.iaSizeY,'%d'),...
                ', grid step ',num2str(pivPar.iaStepX,'%d'), ...
                'x', num2str(pivPar.iaStepY,'%d')]);
            drawnow;
        end
    end
    % remove temporary data from pivData
    pivData = rmfield(pivData,'ccW');
    pivData = orderfields(pivData);
    % get the computational time
    pivData.infCompTime(kp) = toc(timer);
end

% remove temporary data from pivData
pivData = rmfield(pivData,'imArray1');
pivData = rmfield(pivData,'imArray2');
pivData = rmfield(pivData,'imMaskArray1');
pivData = rmfield(pivData,'imMaskArray2');

% add additional fields...
pivData.Nx = size(pivData.X,2);
pivData.Ny = size(pivData.X,1);

% 3. sort pivData's fields alphabeticcally
pivData = orderfields(pivData);

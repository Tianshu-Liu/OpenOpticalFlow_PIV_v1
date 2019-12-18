function [pivPar,pivData] = pivParams(pivData,pivParIn,action,varargin)
% pivParams - adjust content of pivPar variable, which controls setting of PIV analysis
%
% Usage:
%   1. [pivPar,pivData] = pivParams(pivData,pivParIn,'Defaults')
%   2. [pivPar,pivData] = pivParams(pivData,pivParIn,'DefaultsSeq')
%   3. [pivPar,pivData] = pivParams(pivData,pivParIn,'SinglePass',K)
%   4. [pivPar,pivData] = pivParams(pivData,pivParIn,'Defaults1Px')
%   4. [pivPar,pivData] = pivParams(pivData,pivParIn,'DefaultsJobManagement')
%
% Usage 1: Parameters present in pivParIn are copied to output structure pivPar. pivPar is then completed by
% missing parameters, which are set to defaults. This usage is intended for setting pivPar for the use with
% pivAnalyzeImagePair.
%
% Usage 2: Parameters present in pivParIn are copied to output structure pivPar. pivPar is then completed by
% missing parameters, which are set to defaults. This usage is intended for setting pivPar for the use with
% pivAnalyzeImageSequence.
%
% Usage 3: If multipass PIV analysis is performed, pivPar usually contains settings for all passes. Using 
% pivParams(pivData,pivParIn,'SinglePass',K) then separates the parameters for K-th pass.
%
% Usage 4: Parameters present in pivParIn are copied to output structure pivPar. pivPar is then completed by
% missing parameters, which are set to defaults. This usage is intended for setting pivPar for the use with
% pivAnalyzeSequence1Px.
%
% Usage 4: Parameters present in pivParIn are copied to output structure pivPar. pivPar is then completed by
% missing parameters jmParallelJobs and jmLockfileExpirationTime, which are set to defaults. This usage is called
% by pivManageJobs; users usually do not need to use this options themselves.
%
%
% Inputs:
%   pivData ... (struct) Structure with piv data. Can be empty []. If action is 'singlePass', pivPar are added
%          to pivData and send as output.
%   pivParIn ... (struct) Structure with parameters ofr PIV processing. May be empty [] For field names and 
%          meaning see description of pivPar in outputs. 
%   K ... number of PIV pass, for which parameters are required.
%
% Outputs:
%    pivPar ... (struct) parameters defining the evaluation. For action 'Defaults', fields present in 
%          pivParIn will be copied to pivPar. Missing fields will be set to default values. For action
%          'SinglePass', only elements of fields in pivParIn, which are relevant for PassOrder-th pass of PIV,
%          are copied. Fields of pivPar structure are summarized below:
%      --- prefix an - these fields are used mostly in subroutines pivAnalyseImagePair.m
%          anNpasses ... number of passes
%                        - default value: if iaSizeX or iaSizeY is specified, then its length; 4 otherwise
%          anVelocityEst ... way, how velocity estimate is determined before processing each image pair.
%               Parameter affects only analysis of image sequence by pivAnalyzeImageSequence.m. Possible
%               values are
%               'previous' ... velocity field from the previous image pair is taken as the velocity estimate
%               'previousSmooth' ... as previous, but the velocity field is smoothed first
%               'none' ... no velocity field (U = V = 0) is considered
%               'pivData' ... veloccity estimate is taken from pivData structure as the value corresponding
%                    for a given time slice.
%             - default value: 'previous'
%          anOnDrive ... (logical) If true, results of processing of each image pair in a sequence are stored 
%               in the folder specified by .anTargetPath. Affects only pivAnalyzeImageSequence.m.
%             - default value: false
%          anTargetPath ... specified folder, to which files containing results of PIV analysis of each image
%               pair are stored. Affects only pivAnalyzeImageSequence.m.
%             - default value: ''
%          anForceProcessing ... if false and if .anOnDrive is true, processing of an image pair is skipped if
%               an output file for a given image pair already exists. Results stored in the output file is
%               read instead. If .anForceProcessing is true, image pair is analyzed and if the output file
%               exists, it is overwritten. Affects only pivAnalyzeImageSequence.m.
%             - default value: false
%          anPairsOnly ... if true and if .anOnDrive is true, the subroutine pivAnalyzeImageSequence.m 
%               processes only image pairs, but it does not produce final file with sequence data. If this 
%               option is set, much less memory is required for processing.
%          anStatsOnly ... if true, the subroutine pivAnalyzeImageSequence.m processes image pairs and add
%               them to the velocity statistics. The final data of velocity processing includes only velocity
%               statistics, but not complete velocity record. Much less memory is required for processing.
%      --- prefix ia - these fields are used mostly in subroutine pivInterrogate.m
%          iaSizeX, iaSizeY ... size of interrogation area [px]
%             - default value: if the other of iaSizeX, iaSizeY exists, set it to the same value
%             - otherwise, set it folowwing sequence [64, 32, 16, 16, 16, 16, ...] with anNpasses terms
%          iaStepX, iaStepY ... step between interrogation areas [px]
%             - default value: floor(iaSizeX/2) or floor(iaSizeY/2)
%          iaMethod ... way, how interrogation area are created. Possible values are
%               'basic' ... interrogatio areas are regularly distribute rectangles
%               'offset' ... interrogation areas are shifted by the estimated displacement
%               'deflinear' ... deformable interrogation areas with linear deformation
%               'defspline' ... deformable interrogation areas with spline deformation
%             - Note: if Uest and Vest contains only zeros or if they are empty/unspecified, 'basic' method is
%                    always invoked regardless .iaMethod setting
%             - default value: 'defspline'
%          iaImageToDeform ... defines, which image should deform. It is taken in account if .iaMethod is
%                  'deflinear' or 'defspline', or if it is 'offset' (then it defines, in which image IAs are
%                  shifted). Possible values are
%              'image1', 'image2' ... either im1 or im2 is deformed correspondingly to Uest and Vest
%              'both' ... deformation both images are deformed by Uest/2 and Vest/2. Sligthly more CPU time is
%                         required.
%             - default value: 'defspline'
%          iaImageInterpolationMethod ... way, how the images are interpolated when deformable IAs are used
%                  (for .iaMethod == 'deflinear' or 'defspline'. Possible values are:
%              'linear', 'spline' ... interpolation is carried out using interp2 function with option either
%                                     '*linear' or '*spline'
%             - default value: 'natural'
%          iaPreprocMethod ... defines image preprocessing method. Possible values are
%               'none' ... no image preprocessing
%               'MinMax' ... MinMax filter is applied (see p. 248 in Ref. [1])
%             - default value: 'none'
%          iaMinMaxSize ... (applies only if iaPreprocMethod is 'MinMax'). Size of MinMax filter kernel.
%             - default value: 7
%          iaMinMaxLevel ... (applies only if iaPreprocMethod is 'MinMax'). Contrast level, below which
%               contrast in not more enhanced.
%             - default value: 16
%      --- prefix im - relates to the images
%          imMask1, imMask2 ... Masking images for im1 and im2. It should be either empty (no mask), or of the
%               same size as im1 and im2. Masked pixels should be 0 in .imMaskX, non-masked pixels should be 1
%             - default value: empty []
%      --- prefix jm - these fields are related to job management, used mostly by subroutines pivManageJobs.m and
%               pivAnalyzeImageSequemce.m
%          jmParallelJobs ... defines, to how many jobs is the treatment of present task distributed. Used by
%               pivManageJobs.m.
%          jmLockFile ... lock file name and path. If this field exists, this file is rewritten
%               with a status message everytimes an image pair is treated. 
%          jmLockExpirationTime ... maximum age of a lock file. If a lock file is older than this limit, it is
%               regarded as non-existent. Used by pivManageJobs.m.
%      --- prefix cc - these fields are used mostly in subroutines pivCrossCorr.m
%          ccRemoveIAMean ... if =0, do not remove IA's mean before cross-correlation; if =1, remove the mean;
%              if in between, remove the mean partially
%            - default value: 1
%          ccMaxDisplacement ... maximum allowed displacement to accept cross-correlation peak. This parameter
%              is a multiplier; e.g. setting ccMaxDisplacement = 0.6 means that the cross-correlation peak
%              must be located within [-0.6*iaSizeX...0.6*iaSizeX, -0.6*iaSizeY...0.6*iaSizeY] from the zero
%              displacement.
%            - Note: IA offset is not included in the displacement accounted by ccMaxDisplacement, hence real
%              displacement can be larger if ccIAmethod is any other than 'basic'.
%            - default value: 0.9
%          ccWindow ... windowing function. Possible values are 
%              'uniform', 'Gauss', 'parzen', 'Hanning', 'Welch' ...windows described on p. 389 in Ref. [1]
%              'Nogueira'... see p. 396 in Ref. [1]
%              'Hanning2', 'Hanning4' ... Hanning window, which is in second or fourth power, respectively
%              power, respectively (the window is then more narrow).
%            - default value: 'Welch'. If missing in pivPar during evaluation, 'uniform' is considered for
%              compatibility reasons.
%          ccCorrectWindowBias ... Set whether cc is corrected for bias due to IA windowing.
%            - default value: false 
%            - ATTENTION: option true increases strongly noise on data, probably something is wrong.
%          ccMethod ... set methods for finding cross-correlation of interrogation areas. Possible values are
%              'fft' (use Fast Fourrier Transform' and 'dcn' (use discrete convolution). Option 'fft' is more
%              suitable for initial guess of velocity. Option 'dcn' is suitable for final iterations of the
%              velocity field, if the displacement corrections are already small.
%            - default value: 'fft'
%          ccMaxDCNdist ... Defines maximum displacement, for which the correlation is computed by DCN method 
%              (apllies only for ccMethod = 'dcn'). 
%            - default value: 1
%      --- prefix vl - these fields relates to vector validation, subroutine pivValidate.m. Parameters without 
%              "Seq" are used for validation of of individual image pairs. Parameters with "Seq" are used when
%              an image sequence is validated.
%              individual image pairs.
%          vlMinCC ... minimum value of cross-correlation peak. Vectors with .ccPeak < vlMinCC*median(ccPeak) 
%              will be marked as invalid. If vlMinCC == 0, this test is skipped.
%            - default value: 0.3
%          vlTresh, vlEps ... Define treshold for the median test. To accepted, the difference of actual vector
%              from the median (of vectors in the neighborhood) should be at most vlTresh *(vlEps +
%              (neighborhood vectors) - (their median))
%            - default values: vlTresh = 2, vlEps = 0.1
%          vlDist ... to what distance median test is performed (if vlDist = 1, kernel has size 3x3; for
%              vlDist = 2, kernel is 5x5, and so on)
%            - default value: 2
%          vlPasses ... number of passes of the median test
%            - default value: 2
%          vlDistTSeq ... to what distance (in time) the median test is performed. If vlDistT == 1, median test
%              will be based on one previous and one subsequent time slices. This parameter is taken in
%              account only if pivValidate.m is applied on pivData, in which contains data for image sequence.
%            - Default value: 1
%          vlTreshSeq, vlEpsSeq, vlDistSeq,vlPassesSeq ... same as without "Seq" at the end, but used when
%              pivData of an image sequence are validated. 
%            - Default values: vlTreshSeq = 2, vlEpsSeq = 0.06, vlDistSeq = 2, vlPassesSeq = 1
%      --- prefix rp - affects replacement of invalid vectors, subroutine pivReplace.m
%          rpMethod ... specifies how the spurious vectors are replaced. 
%              Possible values are
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
%            - default value: 'inpaintT'
%      --- prefix cr - controls, how predictor-corrector correction is applied (affects pivCorrector.m). 
%          crAmount ... Controls amount, to which predictor-corrector is applied. If 0, correction is not
%              applied (pivCorrector.m quits withoud modifying the velocity field). If 1, correction is
%              applied. If 0 < crAmount < 1, correction is applied, but it is weaker than it should be.
%            - default value: 0
%      --- prefix seq - controls, which images are processed as image pairs when dealing with sequences
%          (affects pivCreateImageaSequence.m). Parameters are:
%          seqPairInterval ... (if unspecified, default value is 2): interval between the index of first 
%              image in consecutive pairs. E.g., if PairInterval == 1, image pairs will be Img01 + Img02, 
%              Img02 + Img03, Img03 + Img04, etc. If PairInterval == 2, image pairs will be Img01 + Img02, 
%              Img03 + Img04, Img05 + Img06, etc. For PairInterval == 5, the sequence is Img01+Img02, 
%              Img06+Img07, Img11+Img12, etc.
%          seqFirstIm ... (if unspecified, default value is 1): position of the first image of the first pair 
%              in the imagelist. E.g., if PairInterval = 2 and FirstIm = 2, image pairs are Img02 + Img03, 
%              Img04 + Img05, Img06 + Img07, etc. 
%          seqDiff ... (if unspecified, default value is 1): difference between index of images within one 
%              pair. Allowed values are 1, 3, 5, 7, ... . E.g. if PairInterval == 5, FirstIm = 5 and Diff = 1, 
%              image pairs are Img05 + Img06, Img10 + Img11, Img15 + Img16. If PairInterval == 5 and Diff = 1, 
%              FirstIm = 5 and Diff = 3, image pairs are Img04 + Img07, Img09 + Img12, Img14 + Img17 etc. This 
%              parameter allows to increase time difference between images in the case of time-resolved 
%              records.
%          seqMaxPairs ... (if unspecified, default value is +Inf): Maximum number of images in the output 
%              lists.
%      --- prefix sm - affects smoothing of vector field, subroutine pivSmooth.m
%          smMethod ... defines smoothing method. Possible values are:
%              'none' ... do not perform smoothing
%              'smoothn' ... uses smoothn.m function by Damian Garcia [5]
%              'gauss' ... uses Gaussian kernel
%            - default value: 'smoothn' for all passes except the last, for which it is 'none'
%          smSigma ... amount of smoothing
%            - default value: 0.2
%          smSize ... size of filter (applies only to Gaussian filter)
%            - default value: 5 (this default value is set only if smMethod = 'gauss')
%          smMethodSeq ... defines smoothing method for processing image sequences. Possible values are:
%              'none' ... do not perform smoothing
%              'smoothn' ... uses smoothn.m function by Damian Garcia [5]
%            - default value: 'none'
%          smSigmaSeq ... amount of smoothing
%            - default value: 1
%      --- prefix sp - affects analysis carried out using single-pixel correlation, subroutine 
%              pivAnalyzeSequence1Px.m and its subroutines
%          spDeltaXNeq, spDeltaXPos ... defines maximum displacement in X direction, for which cross-
%              correlation function is evaluated
%            - default value: 8
%          spDeltaYNeq, spDeltaYPos ... defines maximum displacement in Y direction, for which cross-
%              correlation function is evaluated
%            - default value: 8
%          spDeltaAutoCorr ... maximum displacement in both X and Y directions, for which auto-correlation 
%              function is evaluated
%            - default value: 3
%          spBindX, spBindY ... defines binding of image pixels. For example, if spBindX=2 and spBindY=4, 
%              cross-correlations of 2x4 neighbouring pixels are averaged before velocity field is computed. 
%              Binding decreases resulting resolution.
%            - default value: 1 (no binding)
%          spACsource ... defines, which image is used for computing the autocorrelation function 
%              (pivSinglepixCorrelate.m). Possible values are 
%              'both' ... the autocorralation is evaluated for both first and second image in each pair and 
%                  the average of these two autocorrelations is taken in account; use this option when first 
%                  and second laser pulses are not the same
%              'im1', 'im2' ... only first or second frame in each pair is used for evaluation of the
%                  auto-correlation function.
%            - default value: 'both'
%          spAvgSmooth, spRmsSmooth ... defines level of smoothing of average and rms images, respectively. 
%            - default value: spAvgSmooth = 3, spRmsSmooth = 5
%          spGFitNPasses ... number of passes for fitting of gaussian distribution to the cross-correlation peaks.
%              Results from the previous pass are used (after smoothing) as the initial guess for the optimization.
%            - default value: 1
%          spGFitMinCc ... when fitting cross-correlation peak, displacements, for which cross-correlation is 
%              smaller than spMinFitCC*(max(CC)), are not considered. Affects pivSinglepixGaussFit.m. spGFitMinCc should
%              be 1D array with number of elements given by spGFitNPasses.
%            - default value: 0.2
%          spGFitMaxDist ... defines the maximum distance from the cross-correlation peak, for which the
%              evaluated CC if fitted by the gaussian distribution. Affects pivSinglepixGaussFit.m. spGFitMaxDist should
%              be 1D array with number of elements given by spGFitNPasses.
%            - default value: 5
%          spGFitMaxIter ... determines maximum number of iterations during fitting of gaussian peak.
%            - default value: 10000
%          spVlTresh, spVlEps ... Define treshold for the validation by the median test. To accepted, the 
%              difference of actual values from the median (in the neighborhood) should be at most 
%              vlTresh *(vlEps + (neighborhood values) - (their median)). Affects pivSinglepixValidate.m. These
%              parameters should be 1D array with number of elements given by spGFitNPasses.
%            - default values: spVlTresh = 1.8, spVlEps = 0.05
%          spVlDist ... determines, to what distance the median test is performed (if vlDist = 1, kernel has 
%              size 3x3; for vlDist = 2, kernel is 5x5, and so on). Affects pivSinglepixValidate.m. This
%              parameters should be 1D array with number of elements given by spGFitNPasses.
%            - default value: 2
%          spVlPasses ... number of passes of the median test. Affects pivSinglepixValidate.m. This parameters should be 
%              1D array with number of elements given by spGFitNPasses.
%            - default value: 1
%          spSmSigma ... amount of smoothing of the first fit of CC peak (used as a guess for the second fit). Affects 
%              pivSinglepixSmooth.m. This parameters should be 1D array with number of elements given by spGFitNPasses.
%            - default value: 2 
%          spOnDrive ... determines, whether the processing results are stored on disk.
%            - default value: true (recommended)
%          spForceProcessing ... determines, whether the processing should be performed, even if an output 
%              file with intermediate result exists
%            - default value: false
%          spSaveInterval ... defines, how many images are correlated before storing an intermediate result
%              (pivSinglepixCorrelate.m)
%            - default value: 50
%      --- Following fields (with qv prefix) are used to plot results during computing.
%          qvPair ... (cell array) string of options that is send to pivQuiver.m to plot intermediate
%              results after each processing pass by pivAnalyseImagePair.m
%            - default value: {'Umag','quiver','selectStat','valid','linespec','-k','quiver','selectStat',
%                  'replaced','linespec','-w'} (valid black + white replaced vectors, on background with 
%                  velocity magnitude);
%            - for options for pivQuiver, see help of pivQuiver
%    pivData ... Structure with PIV data. If action is 'SinglePass', parameters pivPar are added to pivData at
%         the input and copied to output
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
%    Seq ... concerns working with a sequence of velocity fields
%    sm ... smoothing
%    sp ... single pixel
%    Word "velocity" should be understood as "displacement"
%%

switch lower(action)
    case 'defaults'  % SET DEFAULT VALUES TO UNSPECIFIED FIELDS OF PARIN
        pivPar = pivParIn;
        % anNpasses - number of passes:
        if ~isfield(pivPar,'anNpasses')
            if isfield(pivPar,'iaSizeX')
                pivPar.anNpasses = numel(pivPar.iaSizeX);
            elseif isfield(pivPar,'iaSizeY')
                pivPar.anNpasses = numel(pivPar.iaSizeY);
            else
                pivPar.anNpasses = 4;
            end
        end
        % iaSizeX, iaSizeY - size of interrogation area:
        if ~isfield(pivPar,'iaSizeX')
            if isfield(pivPar,'iaSizeY')
                pivPar.iaSizeX = pivPar.iaSizeY;
            else
                aux = [64,32];
                if pivPar.anNpasses > numel(aux), aux = [aux, 16*ones(1,pivPar.anNpasses-numel(aux))]; end
                pivPar.iaSizeX = aux(1:pivPar.anNpasses);
            end
        end
        if ~isfield(pivPar,'iaSizeY')
            if isfield(pivPar,'iaSizeX')
                pivPar.iaSizeY = pivPar.iaSizeX;
            else
                aux = [64,32,32];
                if pivPar.anNpasses > 3, aux = [aux, 16*ones(1,pivPar.anNpasses-3)]; end
                pivPar.iaSizeY = aux(1:pivPar.anNpasses);
            end
        end
        if ~isfield(pivPar,'iaStepX')
            if isfield(pivPar,'iaStepY')
                pivPar.iaStepX = pivPar.iaStepY;
            else
                pivPar.iaStepX = floor(pivPar.iaSizeX/2);
            end
        end        % other pivParameters
        if ~isfield(pivPar,'iaStepY')
            if isfield(pivPar,'iaStepX')
                pivPar.iaStepY = pivPar.iaStepX;
            else
                pivPar.iaStepY = floor(pivPar.iaSizeY/2);
            end
        end        % other pivParameters
        pivPar = chkfield(pivPar,'iaMethod','defspline');
        pivPar = chkfield(pivPar,'imMask1',[]);
        pivPar = chkfield(pivPar,'imMask2',[]);
        pivPar = chkfield(pivPar,'iaImageToDeform','image1');
        pivPar = chkfield(pivPar,'iaImageInterpolationMethod','spline');
        pivPar = chkfield(pivPar,'iaPreprocMethod','none');
        pivPar = chkfield(pivPar,'ccRemoveIAMean',1);
        pivPar = chkfield(pivPar,'ccMaxDisplacement',0.9);
        pivPar = chkfield(pivPar,'ccWindow','Welch');
        pivPar = chkfield(pivPar,'ccCorrectWindowBias',false);
        pivPar = chkfield(pivPar,'ccMaxDCNdist',1);
        pivPar = chkfield(pivPar,'crAmount',0);
        pivPar = chkfield(pivPar,'vlMinCC',0.3);
        pivPar = chkfield(pivPar,'vlTresh',2);
        pivPar = chkfield(pivPar,'vlEps',0.1);
        pivPar = chkfield(pivPar,'vlPasses',[2 1 1]);
        pivPar = chkfield(pivPar,'vlDist', 2*ones(pivPar.anNpasses,1));
        pivPar = chkfield(pivPar,'smMethod','smoothn');
        pivPar = chkfield(pivPar,'smSigma',NaN);
        pivPar = chkfield(pivPar,'rpMethod','inpaint');
        pivPar = chkfield(pivPar,'qvPair',{}); 
           % possible setting:{'Umag','quiver','selectStat','valid','linespec','-k','quiver','selectStat',...
           % 'replaced','linespec','-w'}
        % set smSize only if smMethod = Gauss
        aux = false;
        if iscell(pivPar.smMethod)
            for kk = 1:numel(pivPar.smMethod), aux = aux + strcmpi(pivPar.smMethod{kk},'gauss'); end
        elseif ischar(pivPar.smMethod)
            aux = strcmpi(pivPar.smMethod,'gauss');
        end
        if aux
            pivPar = chkfield(pivPar,'smSize',5);
        end
        % set ccMethod to 'dcn' if ia size is 12 or smaller, or if previous iteration used the same IA size
        aux = cell(1,numel(pivPar.iaSizeX));
        for kk = 1:numel(pivPar.iaSizeX)
            aux{kk} = 'dcn';
            if kk==1, aux{kk} = 'fft'; end
            if kk>1 && (max(pivPar.iaSizeX(kk),pivPar.iaSizeY(kk)) > 12) && ...
                ~((pivPar.iaSizeX(kk) == pivPar.iaSizeX(kk-1)) && ...
                    (pivPar.iaSizeY(kk) == pivPar.iaSizeY(kk-1)))
                aux{kk} = 'fft';
            end
        end
        pivPar = chkfield(pivPar,'ccMethod',aux);
        % set validation to have two passes in the first and last PIV pass
        aux = ones(1,pivPar.anNpasses);
        aux(1) = 2; aux(end) = 2;
        pivPar = chkfield(pivPar,'vlPasses',aux);
        % set validation distance to be 1 in the first pass, 2 in other passes
        aux = 2*ones(1,pivPar.anNpasses);
        aux(1) = 1;
        pivPar = chkfield(pivPar,'vlDist',aux);
        % set iaMinMaxSize and iaMinMaxLeve, is applicable
        if strcmpi(pivPar.iaPreprocMethod,'MinMax')
            pivPar = chkfield(pivPar,'iaMinMaxSize',7);
            pivPar = chkfield(pivPar,'iaMinMaxLevel',16);
        end
        % sort fields alphabetically
        pivPar = orderfields(pivPar);
        
        
    case 'defaultsseq'   % SET DEFAULTS WHEN PROCESSING A SEQUENCE OF IMAGES
        % add additional defaults
        pivPar = pivParIn;
        pivPar = chkfield(pivPar,'vlDistTSeq',0);
        pivPar = chkfield(pivPar,'vlTreshSeq',2);
        pivPar = chkfield(pivPar,'vlEpsSeq',0.1);
        pivPar = chkfield(pivPar,'vlPassesSeq',1);
        pivPar = chkfield(pivPar,'vlDistSeq',2);
        pivPar = chkfield(pivPar,'smMethodSeq','none');
        pivPar = chkfield(pivPar,'smSigmaSeq',1);
        pivPar = chkfield(pivPar,'seqMaxPairs',+Inf);
        pivPar = chkfield(pivPar,'seqFirstIm',1);
        pivPar = chkfield(pivPar,'seqDiff',1);
        pivPar = chkfield(pivPar,'seqPairInterval',1);
        pivPar = chkfield(pivPar,'anOnDrive',false);
        pivPar = chkfield(pivPar,'anTargetPath','');
        pivPar = chkfield(pivPar,'anForceProcessing',false);
        pivPar = chkfield(pivPar,'anVelocityEst','previous');
        pivPar = chkfield(pivPar,'anPairsOnly',false);
        pivPar = chkfield(pivPar,'anStatsOnly',false);        
        if strcmpi(pivPar.anVelocityEst,'previous')||strcmpi(pivPar.anVelocityEst,'previousSmooth')
            pivPar = chkfield(pivPar,'ccMethod','dcn');
        end
        % set defaults for each image pair
        [pivPar,pivData] = pivParams(pivData,pivPar,'Defaults');
        % sort fields alphabetically
        pivPar = orderfields(pivPar);

    case 'defaults1px'   % SET DEFAULTS WHEN PROCESSING A SEQUENCE OF IMAGES USING SINGLE-PIXEL CORRELATION
        % set defaults for each image pair
        [pivPar,pivData] = pivParams(pivData,pivParIn,'Defaults');
        % add additional defaults
        pivPar = chkfield(pivPar,'spDeltaXNeg',8);
        pivPar = chkfield(pivPar,'spDeltaXPos',8);
        pivPar = chkfield(pivPar,'spDeltaYNeg',8);
        pivPar = chkfield(pivPar,'spDeltaYPos',8);
        pivPar = chkfield(pivPar,'spDeltaAutoCorr',3);
        pivPar = chkfield(pivPar,'spBindX',1);
        pivPar = chkfield(pivPar,'spBindY',1);
        pivPar = chkfield(pivPar,'spStepX',min(pivPar.spBindX,pivPar.spBindY));
        pivPar = chkfield(pivPar,'spStepY',min(pivPar.spBindX,pivPar.spBindY));
        pivPar = chkfield(pivPar,'spAvgSmooth',3);
        pivPar = chkfield(pivPar,'spRmsSmooth',5);
        pivPar = chkfield(pivPar,'spACsource','both');
        pivPar = chkfield(pivPar,'spOnDrive',true);
        pivPar = chkfield(pivPar,'spForceProcessing',false);
        pivPar = chkfield(pivPar,'spSaveInterval',50);
        pivPar = chkfield(pivPar,'spGFitNPasses',1);
        pivPar = chkfield(pivPar,'spGFitMaxIter',10000);
        if ~isfield(pivPar,'spGFitMinCc')
            pivPar.spGFitMinCc = 0.05*ones(1,pivPar.spGFitNPasses);
        elseif numel(pivPar.spGFitMinCc)==1
            pivPar.spGFitMinCc = pivPar.spGFitMinCc*ones(1,pivPar.spGFitNPasses);
        end
        if ~isfield(pivPar,'spGFitMaxDist')
            pivPar.spGFitMaxDist = 5*ones(1,pivPar.spGFitNPasses);
        elseif numel(pivPar.spGFitMaxDist)==1
            pivPar.spGFitMaxDist = pivPar.spGFitMaxDist*ones(1,pivPar.spGFitNPasses);
        end
        if ~isfield(pivPar,'spVlDist')
            pivPar.spVlDist = 2*ones(1,pivPar.spGFitNPasses);
        elseif numel(pivPar.spVlDist)==1
            pivPar.spVlDist = pivPar.spVlDist*ones(1,pivPar.spGFitNPasses);
        end
        if ~isfield(pivPar,'spVlPasses')
            pivPar.spVlPasses = 2*ones(1,pivPar.spGFitNPasses);
        elseif numel(pivPar.spVlPasses)==1
            pivPar.spVlPasses = pivPar.spVlPasses*ones(1,pivPar.spGFitNPasses);
        end
        if ~isfield(pivPar,'spVlTresh')
            pivPar.spVlTresh = 1.8*ones(1,pivPar.spGFitNPasses);
        elseif numel(pivPar.spVlTresh)==1
            pivPar.spVlTresh = pivPar.spVlTresh*ones(1,pivPar.spGFitNPasses);
        end
        if ~isfield(pivPar,'spVlEps')
            pivPar.spVlEps = 0.05*ones(1,pivPar.spGFitNPasses);
        elseif numel(pivPar.spVlEps)==1
            pivPar.spVlEps = pivPar.spVlEps*ones(1,pivPar.spGFitNPasses);
        end
        if ~isfield(pivPar,'spSmSigma')
            pivPar.spSmSigma = 2*ones(1,pivPar.spGFitNPasses);
            pivPar.spSmSigma(end) = NaN;
        elseif numel(pivPar.spSmSigma)==1
            pivPar.spSmSigma = pivPar.spSmSigma*ones(1,pivPar.spGFitNPasses);
            pivPar.spSmSigma(end) = NaN;
        end
        pivPar = chkfield(pivPar,'seqFirstIm',1);
        pivPar = chkfield(pivPar,'seqDiff',1);
        pivPar = chkfield(pivPar,'seqPairInterval',1);        % sort fields alphabetically
        fieldNames = sort(fieldnames(pivPar));
        aux = [];
        for kk = 1:numel(fieldNames)
            aux.(fieldNames{kk}) = pivPar.(fieldNames{kk});
        end
        pivPar = aux;

    case 'singlepass'   % EXTRACT PAR SPECIFIC TO K-th PASS
        % get pass no.
        kpass = varargin{1};
        % build pivPar for kth pass
        pivPar = [];
        pivPar = copyvalue(pivPar,pivParIn,'anNpasses',kpass,0);
        pivPar = copyvalue(pivPar,pivParIn,'iaSizeX',kpass,2);
        pivPar = copyvalue(pivPar,pivParIn,'iaSizeY',kpass,2);
        pivPar = copyvalue(pivPar,pivParIn,'iaStepX',kpass,2);
        pivPar = copyvalue(pivPar,pivParIn,'iaStepY',kpass,2);
        pivPar = copyvalue(pivPar,pivParIn,'imMask1',kpass,2);
        pivPar = copyvalue(pivPar,pivParIn,'imMask2',kpass,2);
        pivPar = copyvalue(pivPar,pivParIn,'iaMethod',kpass,2);
        pivPar = copyvalue(pivPar,pivParIn,'iaImageToDeform',kpass,2);
        pivPar = copyvalue(pivPar,pivParIn,'iaImageInterpolationMethod',kpass,2);
        pivPar = copyvalue(pivPar,pivParIn,'iaPreprocMethod',kpass,2);
        pivPar = copyvalue(pivPar,pivParIn,'iaMinMaxSize',kpass,2);
        pivPar = copyvalue(pivPar,pivParIn,'iaMinMaxLevel',kpass,2);
        pivPar = copyvalue(pivPar,pivParIn,'ccRemoveIAMean',kpass,2);
        pivPar = copyvalue(pivPar,pivParIn,'ccMaxDisplacement',kpass,2);
        pivPar = copyvalue(pivPar,pivParIn,'ccWindow',kpass,2); 
        pivPar = copyvalue(pivPar,pivParIn,'ccCorrectWindowBias',kpass,2); 
        pivPar = copyvalue(pivPar,pivParIn,'ccMethod',kpass,2); 
        pivPar = copyvalue(pivPar,pivParIn,'ccMaxDCNdist',kpass,2); 
        pivPar = copyvalue(pivPar,pivParIn,'crAmount',kpass,2); 
        pivPar = copyvalue(pivPar,pivParIn,'vlMinCC',kpass,2);
        pivPar = copyvalue(pivPar,pivParIn,'vlTresh',kpass,2);
        pivPar = copyvalue(pivPar,pivParIn,'vlEps',kpass,2);
        pivPar = copyvalue(pivPar,pivParIn,'vlDist',kpass,2);
        pivPar = copyvalue(pivPar,pivParIn,'vlPasses',kpass,2);
        pivPar = copyvalue(pivPar,pivParIn,'vlDistTSeq',kpass,1);
        pivPar = copyvalue(pivPar,pivParIn,'vlTreshSeq',kpass,1);
        pivPar = copyvalue(pivPar,pivParIn,'vlEpsSeq',kpass,1);
        pivPar = copyvalue(pivPar,pivParIn,'vlDistSeq',kpass,1);
        pivPar = copyvalue(pivPar,pivParIn,'vlPassesSeq',kpass,1);
        pivPar = copyvalue(pivPar,pivParIn,'rpMethod',kpass,2);
        pivPar = copyvalue(pivPar,pivParIn,'smMethod',kpass,2);
        pivPar = copyvalue(pivPar,pivParIn,'smSigma',kpass,2);
        pivPar = copyvalue(pivPar,pivParIn,'smSize',kpass,2);
        pivPar = copyvalue(pivPar,pivParIn,'smMethodSeq',kpass,1);
        pivPar = copyvalue(pivPar,pivParIn,'smSigmaSeq',kpass,1);
        pivPar = copyvalue(pivPar,pivParIn,'qvPair',kpass,1);
        % sort fields alphabetically
        fieldNames = sort(fieldnames(pivPar));
        aux = [];
        for kk = 1:numel(fieldNames)
            aux.(fieldNames{kk}) = pivPar.(fieldNames{kk});
        end
        pivPar = aux;
        % store pivParameters in pivData
        if nargout>1
            pivData.pivPar{kpass} = [];
            pivData.pivPar{kpass} = copyvalue(pivData.pivPar{kpass},pivPar,'iaSizeX',kpass,1);
            pivData.pivPar{kpass} = copyvalue(pivData.pivPar{kpass},pivPar,'iaSizeY',kpass,1);
            pivData.pivPar{kpass} = copyvalue(pivData.pivPar{kpass},pivPar,'iaStepX',kpass,1);
            pivData.pivPar{kpass} = copyvalue(pivData.pivPar{kpass},pivPar,'iaStepY',kpass,1);
            % masks are copied, only if pivPar contains a flag for storing extended amount of information
            if isfield(pivParIn,'infoSaveMasks') && pivParIn.infoSaveMasks
                pivData.pivPar{kpass} = copyvalue(pivData.pivPar{kpass},pivPar,'imMask1',kpass,1);
                pivData.pivPar{kpass} = copyvalue(pivData.pivPar{kpass},pivPar,'imMask2',kpass,1);
            else % or store them, if it is text string (path to file with mask)
                if ischar(pivPar.imMask1)
                    pivData.pivPar{kpass} = copyvalue(pivData.pivPar{kpass},pivPar,'imMask1',kpass,1);
                else
                    pivData.pivPar{kpass} = copyvalue(pivData.pivPar{kpass},pivPar,'imMask1',kpass,0);
                end
                if ischar(pivPar.imMask1)
                    pivData.pivPar{kpass} = copyvalue(pivData.pivPar{kpass},pivPar,'imMask2',kpass,1);
                else
                    pivData.pivPar{kpass} = copyvalue(pivData.pivPar{kpass},pivPar,'imMask2',kpass,0);
                end
            end
            % other pivParameters
            pivData.pivPar{kpass} = copyvalue(pivData.pivPar{kpass},pivPar,'ccRemoveIAMean',kpass,1);
            pivData.pivPar{kpass} = copyvalue(pivData.pivPar{kpass},pivPar,'ccMaxDisplacement',kpass,1);
            pivData.pivPar{kpass} = copyvalue(pivData.pivPar{kpass},pivPar,'iaMethod',kpass,1);
            pivData.pivPar{kpass} = copyvalue(pivData.pivPar{kpass},pivPar,'iaImageToDeform',kpass,1);
            pivData.pivPar{kpass} = copyvalue(pivData.pivPar{kpass},pivPar,'iaImageInterpolationMethod',kpass,1);
            pivData.pivPar{kpass} = copyvalue(pivData.pivPar{kpass},pivPar,'rpMethod',kpass,1);
            pivData.pivPar{kpass} = copyvalue(pivData.pivPar{kpass},pivPar,'smMethod',kpass,1);
            pivData.pivPar{kpass} = copyvalue(pivData.pivPar{kpass},pivPar,'smSigma',kpass,1);
            pivData.pivPar{kpass} = copyvalue(pivData.pivPar{kpass},pivPar,'smSize',kpass,1);
            pivData.pivPar{kpass} = copyvalue(pivData.pivPar{kpass},pivPar,'smMethodSeq',kpass,1);
            pivData.pivPar{kpass} = copyvalue(pivData.pivPar{kpass},pivPar,'smSigmaSeq',kpass,1);
            pivData.pivPar{kpass} = copyvalue(pivData.pivPar{kpass},pivPar,'vlTresh',kpass,1);
            pivData.pivPar{kpass} = copyvalue(pivData.pivPar{kpass},pivPar,'vlEps',kpass,1);
            pivData.pivPar{kpass} = copyvalue(pivData.pivPar{kpass},pivPar,'vlDist',kpass,1);
            pivData.pivPar{kpass} = copyvalue(pivData.pivPar{kpass},pivPar,'vlDistT',kpass,1);
            pivData.pivPar{kpass} = copyvalue(pivData.pivPar{kpass},pivPar,'vlPasses',kpass,1);
            pivData.pivPar{kpass} = copyvalue(pivData.pivPar{kpass},pivPar,'vlTreshSeq',kpass,1);
            pivData.pivPar{kpass} = copyvalue(pivData.pivPar{kpass},pivPar,'vlEpsSeq',kpass,1);
            pivData.pivPar{kpass} = copyvalue(pivData.pivPar{kpass},pivPar,'vlDistSeq',kpass,1);
            pivData.pivPar{kpass} = copyvalue(pivData.pivPar{kpass},pivPar,'vlPassesSeq',kpass,1);
        end

        
    case 'defaultsjobmanagement'   % SET DEFAULTS WHEN DEFINING JOBS
        % add additional defaults
        pivPar = pivParIn;
        pivPar = chkfield(pivPar,'jmParallelJobs',4);
        pivPar = chkfield(pivPar,'jmLockExpirationTime',600);
        pivPar = orderfields(pivPar);


end

pivPar = orderfields(pivPar);  % ordder fields alphabetically
end


%% local functions

function [pivPar] = chkfield(pivPar,fieldName,defaultVal)
% check, if field "fieldName" is present in pivParIn. If not present, create it and set its value to defaultVal
if ~isfield(pivPar,fieldName)
    pivPar.(fieldName) = defaultVal;
end
end

function [pivPar] = copyvalue(pivPar,pivParIn,fieldName,k,mode)
% copy field fieldName from pivParIn to pivPar and adapts it to single-pass analysis. Possible copying modes are:
%   0 - do not copy the field
%   1 - copy entire field, if it exists
%   2 - if pivParIn.fieldName is either array or cell, copy only its kth element. If numel(pivParIn.fieldName)<k,
%       copy the last element. If pivParIn.fieldName does not exist, do nothing.
switch mode
    case 0 % do nothing
    case 1
        if isfield(pivParIn,fieldName)
            pivPar.(fieldName) = pivParIn.(fieldName);
        end
    case 2
        if isfield(pivParIn,fieldName)
            if iscell(pivParIn.(fieldName))&& min(size(pivParIn.(fieldName)))==1
                if k<=numel(pivParIn.(fieldName))
                    pivPar.(fieldName) = pivParIn.(fieldName){k};
                else
                    pivPar.(fieldName) = pivParIn.(fieldName){end};
                end
            elseif isnumeric(pivParIn.(fieldName))&& min(size(pivParIn.(fieldName)))==1
                if k<=numel(pivParIn.(fieldName))
                    pivPar.(fieldName) = pivParIn.(fieldName)(k);
                else
                    pivPar.(fieldName) = pivParIn.(fieldName)(end);
                end
            elseif isfield(pivParIn,fieldName)
                pivPar.(fieldName) = pivParIn.(fieldName);
            end
        end
end
end

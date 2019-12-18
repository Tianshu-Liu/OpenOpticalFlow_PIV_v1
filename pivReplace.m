function [pivData] = pivReplace(pivData,pivPar)
% pivReplace - replace displacement vectors, which contains NaN, by values coherent with their neighborhood
%
% Usage:
%     [pivData] = pivReplace(pivData,pivPar)
%
% Inputs:
%     pivData ... (struct) structure containing more detailed results. Required field is
%         X, Y ... position, at which velocity/displacement is calculated
%         U, V ... displacements in x and y direction
%         Status ... matrix describing status of velocity vectors (for values, see Outputs section)
%     pivPar ... (struct) parameters defining the evaluation. Following fields are considered:
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
%
% Outputs:
%    pivData  ... (struct) structure containing more detailed results. If pivData was non-empty at the input, its
%              fields are preserved. Following fields are added or updated:
%        X, Y ... x and y coordinates, at which velocities are evaluated
%        U, V ... x and y components of the velocity/displacement vector (with replaced NaN's)
%        Status ... matrix with statuis of velocity vectors (uint8). Bits have this coding:
%            1 ... masked (set by pivInterrogate)
%            2 ... cross-correlation failed (set by pivCrossCorr)
%            4 ... peak detection failed (set by pivCrossCorr)
%            8 ... indicated as spurious by median test based on image pair (set by pivValidate)
%           16 ... interpolated within one image pair (set by pivReplaced)
%           32 ... smoothed (set by pivSmooth)
%           64 ... indicated as spurious by median test based on image sequence (set by pivValidate)
%          128 ... interpolated within image sequence (set by pivReplaced)
%        replacedN ... number of replaced vectors
%        replacedX,replacedY ... positions, at which velocity/displacement vectors were replaced
%        replacedU,replacedV ... components of the velocity/displacement vectors, which were replaced
%        validN ... number of original and vectors
%        validX,validY ... positions, at which velocity/displacement vectors is original and valid
%        validU,validV ... original and valid components of the velocity/displacement vector
%
%
% This subroutine is a part of
%
% =========================================
%               PIVsuite
% =========================================
%
% PIVsuite is a set of subroutines intended for processing of data acquired with PIV (particle image
% velocimetry).
%
% Written by Jiri Vejrazka, Institute of Chemical Process Fundamentals, Prague, Czech Republic
%
% For the use, see files example_XX_xxxxxx.m, which acompany this file. PIVsuite was tested with
% Matlab 7.12 (R2011a) and 7.14 (R2012a).
%
% In the case of a bug, contact me: vejrazka (at) icpf (dot) cas (dot) cz
%
%
% Requirements:
%     Image Processing Toolbox
%
%     inpaint_nans.m
%         subroutine by John D'Errico, available at http://www.mathworks.com/matlabcentral/fileexchange/4551
%
%     smoothn.m
%         subroutine by Damien Garcia, available at
%         http://www.mathworks.com/matlabcentral/fileexchange/274-smooth
%
% Credits:
%    PIVsuite is a redesigned version of PIVlab software [3], developped by W. Thielicke and E. J. Stamhuis.
%    Some parts of this code are copied or adapted from it (especially from its piv_FFTmulti.m subroutine).
%    PIVsuite uses 3rd party software:
%        inpaint_nans.m, by J. D'Errico, [2]
%        smoothn.m, by Damien Garcia, [5]
%
% References:
%   [1] Adrian & Whesterweel, Particle Image Velocimetry, Cambridge University Press 2011
%   [2] John D'Errico, inpaint_nans subroutine, http://www.mathworks.com/matlabcentral/fileexchange/4551
%   [3] W. Thielicke and E. J. Stamhuid, PIVlab 1.31, http://pivlab.blogspot.com
%   [4] Raffel, Willert, Wereley & Kompenhans, Particle Image Velocimetry: A Practical Guide. 2nd edition,
%       Springer 2007
%   [5] Damien Garcia, smoothn subroutine, http://www.mathworks.com/matlabcentral/fileexchange/274-smooth


% Acronyms and meaning of variables used in this subroutine:
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



%% Replace NaN's in the velocity field

singleType = isa(pivData.U,'single');

X0 = double(pivData.X);
Y0 = double(pivData.Y);
U = double(pivData.U);
V = double(pivData.V);
status = pivData.Status;

method = pivPar.rpMethod;

% if there is only single time slice and method is based also on other time-slices, switch to method for
% single time slice (remov "t" at the end)
if size(U,3) == 1 && (strcmpi(method,'linearT') || strcmpi(method,'naturalT') ...
        || strcmpi(method,'inpaintT'))
    method = method(1:end-1);
end

% detect all elements, for which replacement is required
containsNaN = logical(isnan(U)+isnan(V));
% detect all masked elements (replacement will be removed for them)
masked = logical(bitget(status,1));

switch lower(method)
    case 'none'
        filledU = U;
        filledV = V;
    case {'linear','natural'}
        % create time information
        X = X0;
        Y = Y0;
        filledU = U + NaN;
        filledV = V + NaN;
        for kt = 1:size(U,3)
            % convert everything to vectors
            flatX = reshape(X,numel(X),1);
            flatY = reshape(Y,numel(X),1);
            flatU = reshape(U(:,:,kt),numel(X),1);
            flatV = reshape(V(:,:,kt),numel(X),1);
            % remove all elements, for which U or V is NaN
            OK = logical(~isnan(U).*~isnan(V));
            flatX = flatX(OK);
            flatY = flatY(OK);
            flatU = flatU(OK);
            flatV = flatV(OK);
            % interpolate for all X and Y
            switch lower(method)
                case 'linear'
                    flatUi =  TriScatteredInterp(flatX,flatY,flatU,'linear');
                    flatVi =  TriScatteredInterp(flatX,flatY,flatV,'linear');
                case 'natural'
                    flatUi =  TriScatteredInterp(flatX,flatY,flatU,'natural');
                    flatVi =  TriScatteredInterp(flatX,flatY,flatV,'natural');
            end
            filledU(:,:,kt) = flatUi(X,Y);
            filledV(:,:,kt) = flatVi(X,Y);
        end
    case {'lineart','naturalt'}
        % create time information
        X = U + NaN;
        Y = X;
        T = X;
        if isfield(pivData,'imPairNo')
            for kt = 1:size(T,3)
                T(:,:,kt) = pivData.imPairNo(kt);
                X(:,:,kt) = X0;
                Y(:,:,kt) = Y0;
            end
        else
            for kt = 1:size(T,3)
                T(:,:,kt) = kt;
                X(:,:,kt) = X0;
                Y(:,:,kt) = Y0;
            end
        end
        % convert everything to vectors
        flatX = reshape(X,numel(X),1);
        flatY = reshape(Y,numel(X),1);
        flatT = reshape(T,numel(X),1);
        flatU = reshape(U,numel(X),1);
        flatV = reshape(V,numel(X),1);
        % remove all elements, for which U or V is NaN
        OK = logical(~isnan(U).*~isnan(V));
        flatX = flatX(OK);
        flatY = flatY(OK);
        flatT = flatT(OK);
        flatU = flatU(OK);
        flatV = flatV(OK);
        % interpolate for all X and Y
        switch lower(method)
            case 'lineart'
                flatUi =  TriScatteredInterp(flatX,flatY,flatT,flatU,'linear');
                flatVi =  TriScatteredInterp(flatX,flatY,flatT,flatV,'linear');
            case 'naturalt'
                flatUi =  TriScatteredInterp(flatX,flatY,flatT,flatU,'natural');
                flatVi =  TriScatteredInterp(flatX,flatY,flatT,flatV,'natural');
        end
        filledU = flatUi(X,Y,T);
        filledV = flatVi(X,Y,T);
    case 'inpaint'
        % use D'Errico's inpaiting subroutine for 2D data for each time slice
        filledU = U + NaN;
        filledV = V + NaN;
        for kt = 1:size(U,3)
            filledU(:,:,kt) = inpaint_nans(double(U(:,:,kt)),4);
            filledV(:,:,kt) = inpaint_nans(double(V(:,:,kt)),4);
        end
    case 'inpaintt'
        % use D'Errico's inpaiting subroutine for 3D data for all time slices
        filledU = inpaint_nans3(U,1);
        filledV = inpaint_nans3(V,1);
    case 'inpaintgarcia'
        % use Garcia's inpaiting subroutine for 2D data for each time slice
        filledU = U + NaN;
        filledV = V + NaN;
        for kt = 1:size(U,3)
            filledU(:,:,kt) = inpaintn(U(:,:,kt));
            filledV(:,:,kt) = inpaintn(V(:,:,kt));
        end
    case 'inpaintgarciat'
        % use Garcia's inpaiting subroutine for 3D data for all time slices
        filledU = inpaintn(U);
        filledV = inpaintn(V);
    otherwise
        disp('Error (pivReplace.m): Unknown replacement method.');
        filledU = U;
        filledV = V;
end

% output results
U = filledU;
V = filledV;

% put back NaNs for masked elements
U(masked) = NaN;
V(masked) = NaN;

% output results
replaced = logical(containsNaN .* (~masked));
if singleType
    pivData.U = single(U);
    pivData.V = single(V);
else
    pivData.U = double(U);
    pivData.V = double(V);
end

if size(pivData.U,3) == 1
    status(replaced) = bitset(status(replaced),5);
    pivData.replacedN = sum(sum(replaced));
else
    status(replaced) = bitset(status(replaced),8);
    for kt=1:size(pivData.U,3)
        pivData.replacedN(kt) = sum(sum(replaced(:,:,kt)));
    end
end
pivData.Status = uint16(status);

end


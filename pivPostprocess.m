function [pivDataOut] = pivPostprocess(action,varargin)
% pivPostprocess - postprocess PIV data
%
% Usage:
%    [pivDataOut] = pivPostprocess(action,varargin)
%
%    Typical usage is:
%    1. pivData = pivPostprocess('vorticity',pivData)
%       pivDataSeq = pivPostprocess('vorticity,pivDataSeq)
%          Computes the vorticity field. The vorticity computatio follows the eqs. (9.8) and (9.11) in Adrian
%          & Westerweel, p. 432
%
% Inputs:
%    action ... defines, how the data will be manipulated. Recognized actions (case-insensitive) are
%         'vorticity'.
%         See "Usage" for details.
%    pivData ... structure containing results of pivAnalyzeImagePair.m. In this structure, velocity fields .U
%         and .V have ny x nx elements (where ny and nx is the number of interrogation areas). 
%    pivDataSeq ... structure containing results of pivAnalyzeImageSequence.m. In this structure, velocity
%         fields .U and .V have ny x nx x nt elements (where ny and nx is the number of interrogation areas,
%         nt is the number of image pairs). .ccPeak and .ccPeakSecondary have the same size. Fields
%         .spuriousU, .spuriousV, .spuriousX, .spuriousY are missing.
%
% Outputs:
%    pivData, pivDataSeq ... see inputs for their meaning.
%
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
% Written by Jiri Vejrazka, Institute of Chemical Process Fundamentals, Prague, Czech Republic,
% with these contributions:
%    Nicolas Begue, ENSIACET, Toulouse (action 'computeKuv')
%
% For the use, see files example_XX_xxxxxx.m, which acompany this file. PIVsuite was tested with
% Matlab 7.12 (R2011a) and 7.14 (R2012a).
%
% In the case of a bug, please, contact me: vejrazka (at) icpf (dot) cas (dot) cz
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
%
%
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



switch(lower(action))

    % ACTION computeVorticity
    case {'vorticity'}
        fprintf('Calculating the vorticity... ');
        tic;
        pivDataC = varargin{1};
        vort = zeros(size(pivDataC.U)) + NaN;
        for ki = 1:size(pivDataC.U,3);
            % following. eq. (9.8) and (9.11) in Adrian & Westerweel, p. 432
            % compute filtered velocity fields
            Ufilt = zeros(size(pivDataC.U,1)+2,size(pivDataC.U,2)+2)+NaN;
            Ufilt(2:end-1,2:end-1) = 0.5*pivDataC.U(:,:,ki);
            Ufilt(2:end-1,1:end-2) = Ufilt(2:end-1,1:end-2)+0.25*pivDataC.U(:,:,ki);
            Ufilt(2:end-1,3:end) = Ufilt(2:end-1,3:end)+0.25*pivDataC.U(:,:,ki);
            Vfilt = zeros(size(pivDataC.U,1)+2,size(pivDataC.U,2)+2)+NaN;
            Vfilt(2:end-1,2:end-1) = 0.5*pivDataC.V(:,:,ki);
            Vfilt(1:end-2,2:end-1) = Vfilt(1:end-2,2:end-1)+0.25*pivDataC.V(:,:,ki);
            Vfilt(3:end,2:end-1) = Vfilt(3:end,2:end-1)+0.25*pivDataC.V(:,:,ki);
            % get grid spacing
            dX = pivDataC.X(1,2)-pivDataC.X(1,1);
            dY = pivDataC.Y(2,1)-pivDataC.Y(1,1);
            % compute vorticity
            vort(:,:,ki) = (Vfilt(2:end-1,3:end)-Vfilt(2:end-1,1:end-2))/2/dX - ...
                (Ufilt(3:end,2:end-1)-Ufilt(1:end-2,2:end-1))/2/dY;
        end
        pivDataC.vorticity = single(vort);
        fprintf('finished in %.2f s.\n', toc);
        pivDataOut = pivDataC;
        
        
    % ACTION not recognized
    otherwise
        fprintf('Error: action %s is not recognized by pivPostprocess.m.\n',action);
end
end


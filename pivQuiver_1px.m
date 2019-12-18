function pivQuiver(pivData,varargin)
% pivQuiver - displays quiver plot with the velocity field from PIV analysis with a background image
%
% Usage:
%    1. pivQuiver(pivData,whatToPlot1,option1,option1Value,option2,option2value,...,whatToPlot2,optionN,...)
%          show data for an image pair
%    2. pivQuiver(pivData,'timeSlice',timeSliceNo,whatToPlot1,option1,option1Value,...)
%          choose a timeslice of an image sequence and show data for it
%    3. pivQuiver(pivData,'XYrange',[xmin,xmax,ymin,ymax],whatToPlot1,option1,option1Value,...)
%          choose a part of spacial extent and show data for it
%    4. pivQuiver(pivData,'timeSlice',timeSliceNo,'XYrange',[xmin,xmax,ymin,ymax],whatToPlot1,option1,option1Value,...)
%    5. pivQuiver(pivData,CellWithOptions);
%
% Inputs:
%    pivData ... structure containing results of PIV analysis
%    whatToPlot ... string containing information, what should be shown. Possible values are:
%        'Umag' ... plot background with velocity magnitude
%        'UmagMean' ... plot background with magnitude of an time-average velocity vector
%        'U','V' ... plot background with velocity components
%        'Umean','Vmean' ... plot background with mean values of velocity components
%        'ccPeak' ... plot background with the amplitude of the cross-correlation peak
%        'ccPeakMean' ... plot background with the time-average of amplitude of the cross-correlation peak
%        'ccPeak2nd' ... plot background with the amplitude of the secondary cross-correlation peak
%        'ccPeak2ndMean' ... plot background with the time-average of amplitude of the secondary 
%             cross-correlation peak
%        'ccDetect' ... plot background with the detectability (ccPeak/ccPeakSecondary).
%        'ccDetectMean' ... plot background with the mean detectability (ccPeak/ccPeakSecondary).
%        'ccStd','ccStd1','ccStd2' ... plots background with RMS values in the IA (either average of both
%             images, image1, or image2)
%        'ccMean','ccMean1','ccMean2' ... plots background with mean pixel values in the IA (either average of both
%             images, image1, or image2)
%        'k' ... plot background with turbulence energy k (those should be previously computed using)
%        'RSuu','RSvv','RSuv' ... plot background with Reynolds stress (horizontal normal, vertical normal, or
%             shear)
%        'vort' ... plot background with vorticity (must be computed by pivManipulateData)
%        'epsLEPIV', 'epsMeanLEPIV' ... instantaneous and time-averaged
%             energy dissipation rate, determined by LEPIV (must be computed by pivManipulateData)
%        'image1','image2' ... show image in the background (works only, if images or paths to them are stored
%             in pivData)
%        'imageSup' ... show superposed images in the background (works only, if images or paths to them are stored
%             in pivData)
%        'quiver' ... show quiver plot with the velocity field
%        'invLoc' ... mark spurious vectors with a cross
%        'spQuiver' ... show quiver plot with the mean velocity field obtained by single-pixel correlation
%        'spU','spV' ... plot background with velocity components obtained by single-pixel crosscorrelation;
%            automatically selects between 'spUfit', 'spUint' or 'spU0' (or equivalents for V)
%        'spUmag' ... plot background with mean velocity magnitude obtained by single-pixel crosscorrelation;
%            automatically selects between 'spUfitmag', 'spUintmag' or 'spU0mag'
%        'spU0','spUint','spUfit' ... plot background with velocity components obtained by single-pixel
%            crosscorrelation, with velocity determined by sub-pixel peak interpolation, integration of CC or
%            by optimization, respectively; horizontal velocity component
%        'spV0','spVint','spVfit' ... as option above, but vertical velocity component
%        'spUfilt','spVfilt','spUfiltMag' ... ############# Doplnit #########################
%        'spU0mag','spUintmag','spUfitmag' ... as options above, but plots velocity magnitude
%        'spC1int','spC2int','spPhiint','spC0int' ... results of characterization of CC peak using its
%            integrals
%        'spC1fit','spC2fit','spPhifit','spC0fit' ... results of characterization of CC peak using its fitting
%            by 2D gaussian
%        'spCCmax' ... maximum value of CC
%        'spACC1int','spACC2int','spACPhiint','spACC0int' ... results of characterization of AC peak using its
%            integrals
%        'spACC1fit','spACC2fit','spACPhifit','spACC0fit' ... results of characterization of AC peak using its
%            fitting by 2D gaussian
%        'spACmax' ... maximum value of AC
%        'spCC' ... plot background with expanded single-pixel cross-correlation function. It is recommended 
%            to use this only on a part of image using 'crop' options.
%        'spAC' ... shows expanded single-pixel auto-correlation function. It is recommended to use it only
%            on a part of image using 'crop' options.
%        'spdUdX','spdUdY','spdVdX','spdVdY' ... derivatives of the velocity field ################ doplnit
%        ############
%    optionX ... string defining, which option will be set. Possible values are:
%        'colorMap' ... colormap for background plot. See help for colormap in Matlab's documentation.
%        'clipLo','clipHi' ... set minimum and maximum value of plotted quantity. If applied to quiver plot,
%            applies to velocity magnitude.
%        'subtractU','subtractV' ... velocity, which is subtracted from the velocity field. Typically, mean
%            velocity can be subtracted in order to visualize unsteady structures superposed on the mean flow.
%        'lineSpec' ... (only for quiver) specify appearance of vectors in the quiver plot. Use standard
%            Matlab specifications as in plot command (e.b. '-k' for black solid line, ':y' for dotted yellow
%            vectors etc).
%        'qScale' ... (only for quiver) define length of vectors. If positive, length of vectors is qScale
%            times displacement. If negative, vector length is (-qScale * autoScale), where autoScale is that
%            one for which 20% of vectors are longer than their distance, and 80% are shorter.
%        'crop' ... limits the range of coordinates, for which the plots are shown. Should be followed by
%            array [xmin,xmax,ymin,ymax]. If X or Y should not be limited, adjust it to -Inf or +Inf. This
%            option should be used prior to "WhatToDraw" arguments, as it applies to drawing commands which
%            follows this option.
%        'selectLo','selectHi' ... (only for quiver) ... only vectors whose length is longer than selectLo,
%            and/or shorter than 'selectHi', will be shown.
%        'selectStat' ... (only for quiver) Only vectors with a specific status will be shown. Possible values
%            for status are
%                'valid' ... only non-replaced vectors will be shown
%                'replaced' ... only replaced vectors will be shown
%                'ccFailed' ... only replaced vectors at positions, where cross-correlation failed, will be
%                    shown
%                'invalid' ... only replaced vectors at positions, in which validation indicated spurious
%                    vectors, will be shown
%        'selectMult' ... (only for quiver) only vectors having a specific multiplicator factors will be shown
%        'vecPos' ... (only for quiver) defines, where arrows should be located. If 0, arrows originate in
%            measurement location (center of interrogation area). If 1, arrows finishin the measurement
%            location. If 0.5 (default), arrows midpoint are in measurement positions.
%        'timeSlice' ... (applies only to data containing velocity sequences) Select, which time slice will be
%            considered thereinafter
%        'title' ... writes a title to the plot (should be followed by a string with the title)
%    optionXValue ... value of the optional parameter
%
% Examples:
%    pivQuiver(pivData,'Umag','quiver') ... show background with velocity magnitude and plot velocity vectors
%    pivQuiver(pivData,'Umag','clipHi',5,...
%        'quiver','selectStat','valid','lineSpec','-k',...
%        'quiver','selectStat','replaced','lineSpec','-w') ... show velocity magnitude (clipping longer
%            displacement than 5 pixels), and the velocity field (valid vectors by black, replaced by white
%            color)
%    pivQuiver(pivData,'UmagMean',...
%        'timeSlice',1,'quiver','lineSpec','-k',...
%        'timeSlice',10,'quiver','lineSpec','-y',) ... show mean velocity magnitude and two instantaneous
%            velocity fields at different times, by black and yellow vector

%#ok<*CTCH>

% if no desktop, skip
if ~usejava('jvm') || ~usejava('desktop') || ~feature('ShowFigureWindows')
    return;
end

% get hold status to keep it
holdStatus = ishold;

% get input; can be variable, but also a cell structure
if numel(varargin) == 1 && iscell(varargin{1})
    inputs = varargin{1};
else
    inputs = varargin;
end
ki = 1;

% initializations
pivSlice = pivData;
pivAllSlices = pivData;
sliceNo = NaN;
xmin = -Inf;
xmax = Inf;
ymin = -Inf;
ymax = Inf;

while numel(inputs)>=ki
    if ischar(inputs{ki})
        command = lower(inputs{ki});
        switch command
            case 'timeslice'
                try
                    sliceNo = inputs{ki+1};
                    ki = ki+1;
                    pivSlice = pivManipulateData('readTimeSlice',pivData,sliceNo);
                    if isinf(xmin)&&isinf(xmax)&&isinf(ymin)&&isinf(ymax)
                        pivSlice = pivManipulateData('limitX',pivSlice,[xmin,xmax]);
                        pivSlice = pivManipulateData('limitY',pivSlice,[ymin,ymax]);
                    end
                catch
                    disp('Error (pivQuiver.m): Error when selecting time slice.');
                    return;
                end
                
            case 'crop'
                try
                    auxRange = inputs{ki+1};
                    xmin = auxRange(1);
                    xmax = auxRange(2);
                    ymin = auxRange(3);
                    ymax = auxRange(4);
                    ki = ki+1;
                    if ~isnan(sliceNo)
                        pivSlice = pivManipulateData('readTimeSlice',pivData,sliceNo); 
                    else
                        pivSlice = pivData;
                    end
                    pivSlice = pivManipulateData('limitX',pivSlice,[xmin,xmax]);
                    pivSlice = pivManipulateData('limitY',pivSlice,[ymin,ymax]);
                    pivAllSlices = pivData;
                    pivAllSlices = pivManipulateData('limitX',pivAllSlices,[xmin,xmax]);
                    pivAllSlices = pivManipulateData('limitY',pivAllSlices,[ymin,ymax]);                    
                catch
                    disp('Error (pivQuiver.m): Error when selecting XY range.');
                    return;
                end
            
            case {'umag','umagmean','spumag','spu0mag','spuintmag','spufitmag','spufiltmag',...
                    'u','umean','spu','spu0','spuint','spufit',...
                    'v','vmean','spv','spv0','spvint','spvfit',...
                    'spufilt','spvfilt'}
                % default options
                options.colormap = jet(256);
                options.clipLo = -Inf;
                options.clipHi = +Inf;
                options.subtractU = 0;
                options.subtractV = 0;
                data = pivSlice;
                switch lower(command)
                    case {'umag','u','v'}
                        [options,ki] = parseOptions(inputs,ki+1,options);
                        if strcmpi(options.subtractU,'mean'), options.subtractU = meannan(data.U); end
                        if strcmpi(options.subtractV,'mean'), options.subtractV = meannan(data.V); end
                        qu = (data.U-options.subtractU);
                        qv = (data.V-options.subtractV);
                        xmin = min(data.X);
                        xmax = max(data.X);
                        ymin = min(data.Y);
                        ymax = max(data.Y);
                    case {'umagmean','umean','vmean'}
                        [options,ki] = parseOptions(inputs,ki+1,options);
                        if strcmpi(options.subtractU,'mean'), options.subtractU = meannan(data.U); end
                        if strcmpi(options.subtractV,'mean'), options.subtractV = meannan(data.V); end
                        if isfield(data,'Umean')
                            qu = (data.Umean-options.subtractU);
                            qv = (data.Vmean-options.subtractV);
                        else
                            qu = (data.U-options.subtractU);
                            qv = (data.V-options.subtractV);
                        end
                        xmin = min(data.X);
                        xmax = max(data.X);
                        ymin = min(data.Y);
                        ymax = max(data.Y);
                    case {'spumag','spu','spv'}
                        options.spMinCC = -Inf;
                        [options,ki] = parseOptions(inputs,ki+1,options);
                        if strcmpi(options.subtractU,'mean'), options.subtractU = meannan(data.U); end
                        if strcmpi(options.subtractV,'mean'), options.subtractV = meannan(data.V); end
                        if sum(sum(~isnan(data.spUfit)))>1
                            qu = (data.spUfit-options.subtractU);
                            qv = (data.spVfit-options.subtractV);
                        elseif sum(sum(~isnan(data.spUint)))>1
                            qu = (data.spUint-options.subtractU);
                            qv = (data.spVint-options.subtractV);
                        else
                            qu = (data.spU0-options.subtractU);
                            qv = (data.spV0-options.subtractV);
                        end
                        Ccmax = data.spCCmax;
                        xmin = min(data.spX);
                        xmax = max(data.spX);
                        ymin = min(data.spY);
                        ymax = max(data.spY);
                        auxNOK = logical(Ccmax<options.spMinCC);
                        qu(auxNOK) = NaN;
                        qv(auxNOK) = NaN;
                    case {'spu0mag','spu0','spv0'}
                        options.spMinCC = -Inf;
                        [options,ki] = parseOptions(inputs,ki+1,options);
                        if strcmpi(options.subtractU,'mean'), options.subtractU = meannan(data.U); end
                        if strcmpi(options.subtractV,'mean'), options.subtractV = meannan(data.V); end
                        qu = (data.spU0-options.subtractU);
                        qv = (data.spV0-options.subtractV);
                        Ccmax = data.spCCmax;
                        xmin = min(data.spX);
                        xmax = max(data.spX);
                        ymin = min(data.spY);
                        ymax = max(data.spY);
                        auxNOK = logical(Ccmax<options.spMinCC);
                        qu(auxNOK) = NaN;
                        qv(auxNOK) = NaN;
                    case {'spu0int','spuint','spvint'}
                        options.spMinCC = -Inf;
                        [options,ki] = parseOptions(inputs,ki+1,options);
                        if strcmpi(options.subtractU,'mean'), options.subtractU = meannan(data.U); end
                        if strcmpi(options.subtractV,'mean'), options.subtractV = meannan(data.V); end
                        qu = (data.spUint-options.subtractU);
                        qv = (data.spVint-options.subtractV);
                        Ccmax = data.spCCmax;
                        xmin = min(data.spX);
                        xmax = max(data.spX);
                        ymin = min(data.spY);
                        ymax = max(data.spY);
                        auxNOK = logical(Ccmax<options.spMinCC);
                        qu(auxNOK) = NaN;
                        qv(auxNOK) = NaN;
                    case {'spumagfit','spufit','spvfit'}
                        options.spMinCC = -Inf;
                        [options,ki] = parseOptions(inputs,ki+1,options);
                        if strcmpi(options.subtractU,'mean'), options.subtractU = meannan(data.U); end
                        if strcmpi(options.subtractV,'mean'), options.subtractV = meannan(data.V); end
                        qu = (data.spUfit-options.subtractU);
                        qv = (data.spVfit-options.subtractV);
                        Ccmax = data.spCCmax;
                        xmin = min(data.spX);
                        xmax = max(data.spX);
                        ymin = min(data.spY);
                        ymax = max(data.spY);
                        auxNOK = logical(Ccmax<options.spMinCC);
                        qu(auxNOK) = NaN;
                        qv(auxNOK) = NaN;
                    case {'spufiltmag','spufilt','spvfilt'}
                        options.spMinCC = -Inf;
                        [options,ki] = parseOptions(inputs,ki+1,options);
                        if strcmpi(options.subtractU,'mean'), options.subtractU = meannan(data.U); end
                        if strcmpi(options.subtractV,'mean'), options.subtractV = meannan(data.V); end
                        qu = (data.spUfiltered-options.subtractU);
                        qv = (data.spVfiltered-options.subtractV);
                        Ccmax = data.spCCmax;
                        xmin = min(data.spX);
                        xmax = max(data.spX);
                        ymin = min(data.spY);
                        ymax = max(data.spY);
                        auxNOK = logical(Ccmax<options.spMinCC);
                        qu(auxNOK) = NaN;
                        qv(auxNOK) = NaN;                        
                end
                % calculate Umag and clip it
                if size(qu,3)>1 && (strcmpi(command,'umag')||strcmpi(command,'u')||strcmpi(command,'v'))
                    disp('pivQuiver: Warning: data contains multiple time slices. Plotting mean value.');
                    qu = mean(qu,3);
                    qv = mean(qv,3);
                elseif strcmpi(command,'umagmean')||strcmpi(command,'umean')||strcmpi(command,'vmean')
                    qu = mean(qu,3);
                    qv = mean(qv,3);
                end
                switch lower(command)
                    case {'umag','umagmean','spumag','spu0mag','spuintmag','spufitmag','spufiltmag'}
                        q = sqrt(qu.^2+qv.^2);
                    case {'u','spu','umean','spu0','spuint','spufit','spufilt'}
                        q = qu;
                    case {'v','spv','vmean','spv0','spvint','spvfit','spvfilt'}
                        q = qv;
                end
                q(logical(q < options.clipLo)) = options.clipLo;
                q(logical(q > options.clipHi)) = options.clipHi;
                hold off;
                if options.clipLo == -Inf, qMin = min(min(q)); else qMin = options.clipLo; end
                if options.clipHi == +Inf, qMax = max(max(q)); else qMax = options.clipHi; end
                xmin = min(xmin);
                xmax = max(xmax);
                try
                    imagesc([xmin,xmax],[ymin,ymax],q,[qMin,qMax]);
                catch
                    fprintf('pivQuiver: Failed when plotting %s\n',command);
                end
                axis equal;
                colormap(options.colormap);
                colorbar;
                hold on;
                
            case {'ccpeak','ccpeakmean','ccpeak2nd','ccpeak2ndmean','ccdetect','ccdetectmean',...
                    'ccstd','ccstd1','ccstd2','ccmean','ccmean1','ccmean2',...
                    'k','vort','epslepiv','epsmeanlepiv',...
                    'rsuu','rsvv','rsuv',...
                    'spc1int','spc2int','spphiint','spc0int',...
                    'spc1fit','spc2fit','spphifit','spc0fit',...
                    'spccmax',...
                    'spacc1int','spacc2int','spacphiint','spacc0int',...
                    'spacc1fit','spacc2fit','spacphifit','spacc0fit',...
                    'spacmax',...
                    'spd','spp1','spp2','spc1','spc2',...
                    'spdudx','spdudy','spdvdx','spdvdy',...
                    'sprsuu','sprsvv','sprsuv'}
                % default options
                if strcmpi(command,'vort')
                    options.colormap = vorticityColorMap;
                else
                    options.colormap = jet(256);
                end
                options.clipLo = -Inf;
                options.clipHi = +Inf;
                options.spMinCC = -Inf;
                [options,ki] = parseOptions(inputs,ki+1,options);
                % set quantity q, which will be drawn
                try
                    switch lower(command)
                        case 'ccpeak'
                            data = pivSlice;
                            q = data.ccPeak;
                            if size(q,3)>1
                                disp('pivQuiver: Warning: data contains multiple time slices. Plotting mean value.');
                                q = mean(q,3);
                            end
                        case 'ccpeak2nd'
                            data = pivSlice;
                            q = data.ccPeakSecondary;
                            if size(q,3)>1
                                disp('pivQuiver: Warning: data contains multiple time slices. Plotting mean value.');
                                q = mean(q,3);
                            end
                        case 'ccpeakmean'
                            data = pivAllSlices;
                            q = data.ccPeak;
                            q = mean(q,3);
                        case 'ccpeak2ndmean'
                            data = pivAllSlices;
                            q = data.ccPeakSecondary;
                            aux = isnan(q);
                            q(aux) = mean(q(~aux));
                            q = mean(q,3);
                        case 'ccdetect'
                            data = pivSlice;
                            q1 = data.ccPeak;
                            q2 = data.ccPeakSecondary;
                            aux = logical(isnan(q2));
                            q2(aux) = q1(aux);
                            q = q1./q2;
                            if size(q,3)>1
                                disp('pivQuiver: Warning: data contains multiple time slices. Plotting mean value.');
                                q = mean(q,3);
                            end
                        case 'ccdetectmean'
                            data = pivAllSlices;
                            q1 = data.ccPeak;
                            q2 = data.ccPeakSecondary;
                            aux = isnan(q2);
                            q2(aux) = mean(q2(~aux));
                            q = q1./q2;
                            q = mean(q,3);
                        case 'ccstd'
                            data = pivAllSlices;
                            q = sqrt(data.ccStd1 .* data.ccStd2);
                            q = mean(q,3);
                        case 'ccstd1'
                            data = pivAllSlices;
                            q = data.ccStd1;
                            q = mean(q,3);
                        case 'ccstd2'
                            data = pivAllSlices;
                            q = data.ccStd2;
                            q = mean(q,3);
                        case 'ccmean'
                            data = pivAllSlices;
                            q = 1/2*(data.ccMean1 + data.ccMean2);
                            q = mean(q,3);
                        case 'ccmean1'
                            data = pivAllSlices;
                            q = data.ccMean1;
                            q = mean(q,3);
                        case 'ccmean2'
                            data = pivAllSlices;
                            q = data.ccMean2;
                            q = mean(q,3);
                        case 'k'
                            data = pivAllSlices;
                            q = data.k;
                        case 'vort'
                            try
                                data = pivSlice;
                                q = data.vorticity;
                            catch
                                data = pivPostprocess('vorticity',pivSlice);
                                q = data.vorticity;
                            end
                        case 'epslepiv'
                            try
                                data = pivSlice;
                                q = data.epsLEPIV;
                            catch
                                data = pivPostprocess('LEPIVdissip',pivSlice);
                                q = data.epsLEPIV;
                            end
                        case 'epsmeanlepiv'
                            try
                                data = pivSlice;
                                q = data.epsMeanLEPIV;
                            catch
                                data = pivPostprocess('LEPIVdissip',pivSlice);
                                q = data.epsMeanLEPIV;
                            end                            
                        case 'rsuu'
                            data = pivAllSlices;
                            q = data.RSuu;
                        case 'rsvv'
                            data = pivAllSlices;
                            q = data.RSvv;
                        case 'rsuv'
                            data = pivAllSlices;
                            q = data.RSuv;
                        case 'spc1int'
                            data = pivAllSlices;
                            q = data.spC1int;
                            Ccmax = data.spCCmax;
                            auxNOK = logical(Ccmax<options.spMinCC);
                            q(auxNOK) = NaN;
                        case 'spc2int'
                            data = pivAllSlices;
                            q = data.spC2int;
                            Ccmax = data.spCCmax;
                            auxNOK = logical(Ccmax<options.spMinCC);
                            q(auxNOK) = NaN;
                        case 'spphiint'
                            data = pivAllSlices;
                            q = data.spPhiint;
                            Ccmax = data.spCCmax;
                            auxNOK = logical(Ccmax<options.spMinCC);
                            q(auxNOK) = NaN;
                        case 'spc0int'
                            data = pivAllSlices;
                            q = data.spC0int;
                            Ccmax = data.spCCmax;
                            auxNOK = logical(Ccmax<options.spMinCC);
                            q(auxNOK) = NaN;
                        case 'spc1fit'
                            data = pivAllSlices;
                            q = data.spC1fit;
                            Ccmax = data.spCCmax;
                            auxNOK = logical(Ccmax<options.spMinCC);
                            q(auxNOK) = NaN;
                        case 'spc2fit'
                            data = pivAllSlices;
                            q = data.spC2fit;
                            Ccmax = data.spCCmax;
                            auxNOK = logical(Ccmax<options.spMinCC);
                            q(auxNOK) = NaN;
                        case 'spphifit'
                            data = pivAllSlices;
                            q = data.spPhifit;
                            Ccmax = data.spCCmax;
                            auxNOK = logical(Ccmax<options.spMinCC);
                            q(auxNOK) = NaN;
                        case 'spc0fit'
                            data = pivAllSlices;
                            q = data.spC0fit;
                            Ccmax = data.spCCmax;
                            auxNOK = logical(Ccmax<options.spMinCC);
                            q(auxNOK) = NaN;
                        case 'spccmax'
                            data = pivAllSlices;
                            q = data.spCCmax;
                            Ccmax = data.spCCmax;
                            auxNOK = logical(Ccmax<options.spMinCC);
                            q(auxNOK) = NaN;
                        case 'spacc1int'
                            data = pivAllSlices;
                            q = data.spACC1int;
                            Ccmax = data.spACmax;
                            auxNOK = logical(Ccmax<options.spMinCC);
                            q(auxNOK) = NaN;
                        case 'spacc2int'
                            data = pivAllSlices;
                            q = data.spACC2int;
                            Ccmax = data.spACmax;
                            auxNOK = logical(Ccmax<options.spMinCC);
                            q(auxNOK) = NaN;
                        case 'spacphiint'
                            data = pivAllSlices;
                            q = data.spACPhiint;
                            Ccmax = data.spACmax;
                            auxNOK = logical(Ccmax<options.spMinCC);
                            q(auxNOK) = NaN;
                        case 'spacc0int'
                            data = pivAllSlices;
                            q = data.spACC0int;
                            Ccmax = data.spACmax;
                            auxNOK = logical(Ccmax<options.spMinCC);
                            q(auxNOK) = NaN;
                        case 'spacc1fit'
                            data = pivAllSlices;
                            q = data.spACC1fit;
                            Ccmax = data.spACmax;
                            auxNOK = logical(Ccmax<options.spMinCC);
                            q(auxNOK) = NaN;
                        case 'spacc2fit'
                            data = pivAllSlices;
                            q = data.spACC2fit;
                            Ccmax = data.spACmax;
                            auxNOK = logical(Ccmax<options.spMinCC);
                            q(auxNOK) = NaN;
                        case 'spacphifit'
                            data = pivAllSlices;
                            q = data.spACPhifit;
                            Ccmax = data.spACmax;
                            auxNOK = logical(Ccmax<options.spMinCC);
                            q(auxNOK) = NaN;
                        case 'spacc0fit'
                            data = pivAllSlices;
                            q = data.spACC0fit;
                            Ccmax = data.spACmax;
                            auxNOK = logical(Ccmax<options.spMinCC);
                            q(auxNOK) = NaN;
                        case 'spacmax'
                            data = pivAllSlices;
                            q = data.spACmax;
                            Ccmax = data.spACmax;
                            auxNOK = logical(Ccmax<options.spMinCC);
                            q(auxNOK) = NaN;
                        case 'spdudx'
                            data = pivAllSlices;
                            q = data.spdUdX;
                            Ccmax = data.spCCmax;
                            auxNOK = logical(Ccmax<options.spMinCC);
                            q(auxNOK) = NaN;
                        case 'spdudy'
                            data = pivAllSlices;
                            q = data.spdUdY;
                            Ccmax = data.spCCmax;
                            auxNOK = logical(Ccmax<options.spMinCC);
                            q(auxNOK) = NaN;
                        case 'spdvdx'
                            data = pivAllSlices;
                            q = data.spdVdX;
                            Ccmax = data.spCCmax;
                            auxNOK = logical(Ccmax<options.spMinCC);
                            q(auxNOK) = NaN;
                        case 'spdvdy'
                            data = pivAllSlices;
                            q = data.spdVdY;
                            Ccmax = data.spCCmax;
                            auxNOK = logical(Ccmax<options.spMinCC);
                            q(auxNOK) = NaN;
                        case 'spd'
                            data = pivAllSlices;
                            q = data.spD;
                            Ccmax = data.spCCmax;
                            auxNOK = logical(Ccmax<options.spMinCC);
                            q(auxNOK) = NaN;
                        case 'spc1'
                            data = pivAllSlices;
                            q = data.spC1;
                            Ccmax = data.spCCmax;
                            auxNOK = logical(Ccmax<options.spMinCC);
                            q(auxNOK) = NaN;                            
                        case 'spc2'
                            data = pivAllSlices;
                            q = data.spC2;
                            Ccmax = data.spCCmax;
                            auxNOK = logical(Ccmax<options.spMinCC);
                            q(auxNOK) = NaN;
                        case 'spp1'
                            data = pivAllSlices;
                            q = data.spP1;
                            Ccmax = data.spCCmax;
                            auxNOK = logical(Ccmax<options.spMinCC);
                            q(auxNOK) = NaN;
                        case 'spp2'
                            data = pivAllSlices;
                            q = data.spP2;
                            Ccmax = data.spCCmax;
                            auxNOK = logical(Ccmax<options.spMinCC);
                            q(auxNOK) = NaN;
                        case 'sprsuu'
                            data = pivAllSlices;
                            q = data.spRSuu;
                            Ccmax = data.spCCmax;
                            auxNOK = logical(Ccmax<options.spMinCC);
                            q(auxNOK) = NaN;
                        case 'sprsvv'
                            data = pivAllSlices;
                            q = data.spRSvv;
                            Ccmax = data.spCCmax;
                            auxNOK = logical(Ccmax<options.spMinCC);
                            q(auxNOK) = NaN;
                        case 'sprsuv'
                            data = pivAllSlices;
                            q = data.spRSuv;
                            Ccmax = data.spCCmax;
                            auxNOK = logical(Ccmax<options.spMinCC);
                            q(auxNOK) = NaN;                       
                    end
                    % get min and max of X and Y
                    switch lower(command)
                        case {'ccpeak','ccpeak2nd','ccpeakmean','ccpeak2ndmean','ccdetect','ccdetectmean',...
                                'ccstd','ccstd1','ccstd2','ccmean','ccmean1','ccmean2','k','vort',...
                                'epslepiv','epsmeanlepiv','rsuu','rsvv','rsuv'}
                            xmin = min(min(data.X));
                            xmax = max(max(data.X));
                            ymin = min(min(data.Y));
                            ymax = max(max(data.Y));
                        case {'spc1int','spc2int','spphiint','spc0int','spc1fit','spc2fit','spphifit','spc0fit',...
                              'spacc1int','spacc2int','spacphiint','spacc0int','spacc1fit','spacc2fit',...
                              'spacphifit','spacc0fit','spccmax','spacmax',...
                              'spdudx','spdudy','spdvdx','spdvdy',...
                              'spd','spc1','spc2','spp1','spp2',...
                              'sprsuu','sprsvv','sprsuv'}
                            xmin = min(min(data.spX));
                            xmax = max(max(data.spX));
                            ymin = min(min(data.spY));
                            ymax = max(max(data.spY));
                    end
                catch
                    fprintf('Error (pivQuiver.m): Failed to plot desired data (%s).\n',command);
                    return
                end
                % remove complex numbers
                auxNOK = ~isreal(q);
                if sum(sum(auxNOK))>0
                    fprintf('WARNING (pivQuiver.m): Complex quantity encountered when plotting %s.\n',command);
                    q = abs(q);
                end
                % clip data
                q(logical(q < options.clipLo)) = options.clipLo;
                q(logical(q > options.clipHi)) = options.clipHi;
                hold off;
                if options.clipLo == -Inf, qMin = min(min(q)); else qMin = options.clipLo; end
                if options.clipHi == +Inf, qMax = max(max(q)); else qMax = options.clipHi; end
                % plot quantity q
                imagesc([xmin,xmax],[ymin,ymax],q,[qMin,qMax]);
                axis equal;
                colormap(options.colormap);
                colorbar;
                hold on;
                
            case {'spcc','spac','spccnorm','spccnorm2'}    
                % default options
                options.colormap = jet(256);
                options.colormap(1,:) = 0;
                options.clipLo = -Inf;
                options.clipHi = +Inf;
                [options,ki] = parseOptions(inputs,ki+1,options);
                data = pivAllSlices;
                switch lower(command)
                    case {'spcc','spccnorm','spccnorm2'}
                        q = data.spCC;
                        dXNeg = data.spDeltaXNeg;
                        dXPos = data.spDeltaXPos;
                        dYNeg = data.spDeltaYNeg;
                        dYPos = data.spDeltaYPos;
                    case 'spac'
                        q = data.spAC;
                        dXNeg = (size(q,3)-1)/2;
                        dXPos = dXNeg;
                        dYNeg = dXNeg;
                        dYPos = dXNeg;
                end
                Qflat = zeros(size(q,1)*(size(q,3)+1)+1,size(q,2)*(size(q,4)+1)+1)+NaN;
                for ky = 1:size(q,1)
                    for kx = 1:size(q,2)
                        switch lower(command)
                            case{'spcc','spac'}
                                Qflat(2+(ky-1)*(size(q,3)+1):(ky)*(size(q,3)+1),2+(kx-1)*(size(q,4)+1):(kx)*(size(q,4)+1)) = ...
                                    squeeze(q(ky,kx,:,:));
                            case 'spccnorm'
                                norm0 = min(min(squeeze(q(ky,kx,:,:))));
                                norm1 = max(max(squeeze(q(ky,kx,:,:))))-min(min(squeeze(q(ky,kx,:,:))));
                                if norm1 < 0.01, norm1 = 0.01; end;
                                Qflat(2+(ky-1)*(size(q,3)+1):(ky)*(size(q,3)+1),2+(kx-1)*(size(q,4)+1):(kx)*(size(q,4)+1)) = ...
                                    (squeeze(q(ky,kx,:,:))-norm0)/norm1;
                            case 'spccnorm2'
                                norm1 = max(max(squeeze(q(ky,kx,:,:))))-min(min(squeeze(q(ky,kx,:,:))));
                                if norm1 < 0.01, norm1 = 0.01; end;
                                Qflat(2+(ky-1)*(size(q,3)+1):(ky)*(size(q,3)+1),2+(kx-1)*(size(q,4)+1):(kx)*(size(q,4)+1)) = ...
                                    (squeeze(q(ky,kx,:,:)))/norm1;
                        end
                    end
                end
                Qflat(logical(Qflat < options.clipLo)) = options.clipLo;
                Qflat(logical(Qflat > options.clipHi)) = options.clipHi;
                hold off;
                if options.clipLo == -Inf, qMin = min(min(Qflat)); else qMin = options.clipLo; end
                if options.clipHi == +Inf, qMax = max(max(Qflat)); else qMax = options.clipHi; end
                dx = data.spX(2,2)-data.spX(1,1);
                dy = data.spY(2,2)-data.spY(1,1);
                % imagesc([min(data.spX)-dx/2,max(data.spX)+dx/2],[min(data.spY)-dy/2,max(data.spY)+dy/2],Qflat,[qMin-(qMax-qMin)/254,qMax]);
                imagesc([min(min(data.spX))-dx*(dXNeg+1)/(dXNeg+dXPos+2),max(max(data.spX))+dx*(dXPos+1)/(dXNeg+dXPos+2)],...
                    [min(min(data.spY))-dy*(dYNeg+1)/(dYNeg+dYPos+2),max(max(data.spY))+dy*(dYPos+1)/(dYNeg+dYPos+2)],...
                    Qflat,...
                    [qMin-(qMax-qMin)/254,qMax]);
                axis equal;
                colormap(options.colormap);
                colorbar;
                hold on;

            case {'sppeak0','sppeakint','sppeakfit','spzero'}
                options.lineSpec = '-w';
                options.spMinCC = -Inf;
                data = pivAllSlices;
                N = numel(data.spX);
                X = reshape(data.spX,N,1);
                Y = reshape(data.spY,N,1);
                ccMax = reshape(data.spCCmax,N,1);
                dXNeg = data.spDeltaXNeg;
                dXPos = data.spDeltaXPos;
                dYNeg = data.spDeltaYNeg;
                dYPos = data.spDeltaYPos;
                dx = data.spX(2,2)-data.spX(1,1);
                dy = data.spY(2,2)-data.spY(1,1);
                switch lower(command)
                    case 'spzero'
                        options.pointSpec = 'xy';
                        [options,ki] = parseOptions(inputs,ki+1,options);
                        U = zeros(N,1);
                        V = zeros(N,1);
                    case 'sppeak0'
                        options.pointSpec = 'm.';
                        [options,ki] = parseOptions(inputs,ki+1,options);
                        U = reshape(data.spU0,N,1);
                        V = reshape(data.spV0,N,1);
                    case 'sppeakint'
                        options.pointSpec = 'om';
                        options.lineParts = 60;
                        [options,ki] = parseOptions(inputs,ki+1,options);
                        U = reshape(data.spUint,N,1);
                        V = reshape(data.spVint,N,1);
                        Uel = reshape(data.spU0,N,1);
                        Vel = reshape(data.spV0,N,1);
                        c1 = reshape(data.spC1int,N,1);
                        c2 = reshape(data.spC2int,N,1);
                        phi = reshape(data.spPhiint,N,1);
                    case 'sppeakfit'
                        options.pointSpec = 'xm';
                        options.lineParts = 60;
                        [options,ki] = parseOptions(inputs,ki+1,options);
                        U = reshape(data.spUfit,N,1);
                        V = reshape(data.spVfit,N,1);
                        Uel = reshape(data.spUfit,N,1);
                        Vel = reshape(data.spVfit,N,1);
                        c1 = reshape(data.spC1fit,N,1);
                        c2 = reshape(data.spC2fit,N,1);
                        phi = reshape(data.spPhifit,N,1);
                end
                auxNOK = logical(ccMax<options.spMinCC);
                CCcenterX = X+U*dx/(dXNeg+dXPos+2);
                CCcenterY = Y+V*dy/(dYNeg+dYPos+2);
                CCcenterX(auxNOK) = NaN;
                CCcenterY(auxNOK) = NaN;
                plot(CCcenterX,CCcenterY,options.pointSpec);
                switch lower(command)
                    case {'sppeakint','sppeakfit'}
                        phiel = (0:2*pi/(options.lineParts-1):2*pi)';
                        Xel = zeros(numel(X)*(numel(phiel)+1)+1,1)+NaN;
                        Yel = Xel;
                        for kk=1:numel(X)
                            ksiel = 0.294*c1(kk)*cos(phiel);
                            etael = 0.294*c2(kk)*sin(phiel);
                            xel = Uel(kk)+ksiel.*cos(phi(kk))+etael.*sin(phi(kk));
                            yel = Vel(kk)-ksiel.*sin(phi(kk))+etael.*cos(phi(kk));
                            auxNOK = (xel<-dXNeg-0.5)|(xel>dXPos+0.5)|(yel<-dYNeg-0.5)|(yel>dYPos+0.5);
                            xel(auxNOK) = NaN;
                            yel(auxNOK) = NaN;
                            xel = X(kk) + xel*dx/(dXNeg+dXPos+2);
                            yel = Y(kk) + yel*dy/(dYNeg+dYPos+2);
                            Xel(1+(kk-1)*(numel(phiel)+1):kk*(numel(phiel)+1)-1)=xel;
                            Yel(1+(kk-1)*(numel(phiel)+1):kk*(numel(phiel)+1)-1)=yel;
                        end
                        plot(Xel,Yel,options.lineSpec);
                end
                
                
            case {'image1','image2','imagesup'}
                % default options
                options.colormap = gray(256);
                options.expScale = 1;
                [options,ki] = parseOptions(inputs,ki+1,options);
                try
                    switch lower(command)
                        case 'image1'
                            if ischar(pivSlice.imFilename1)
                                img = imread(pivSlice.imFilename1);
                            elseif isstruct(pivSlice.imFilename1)
                                img = imread(pivSlice.imFilename1{1});
                            end
                        case 'image2'
                            if ischar(pivSlice.imFilename2)
                                img = imread(pivSlice.imFilename2);
                            elseif isstruct(pivSlice.imFilename2)
                                img = imread(pivSlice.imFilename2{1});
                            end
                        case 'imagesup'
                            if ischar(pivSlice.imFilename1)
                                img1 = imread(pivSlice.imFilename1);
                                img2 = imread(pivSlice.imFilename2);
                            elseif isstruct(pivSlice.imFilename1)
                                img1 = imread(pivSlice.imFilename1{1});
                                img2 = imread(pivSlice.imFilename2{1});
                            end
                            img = max(img1,img2);
                    end
                    if isinf(xmin), xmin = 0; end
                    if isinf(xmax), xmax = size(img,2)-1; end
                    if isinf(ymin), ymin = 0; end
                    if isinf(ymax), ymax = size(img,1)-1; end
                    img = img(ymin+1:ymax+1,xmin+1:xmax+1);
                    hold off;
                    image((ymin:ymax)*options.expScale,(xmin:xmax)*options.expScale,img);
                    axis equal;
                    colormap(options.colormap);
                    colorbar('off');
                    hold on;
                catch
                    disp('Warning (pivQuiver.m): Unable to read or display image(s).')
                end
                
            case 'invloc'
                % show by a marker all locations, where i) crosscorelation failed to find a CC peak, ii) where
                % a velocity was marked as spurious
                options.lineSpec = 'xk';
                options.markerSize = 4;
                [options,ki] = parseOptions(inputs,ki+1,options);
                data = pivSlice;
                X = []; Y = [];
                % get list of spurious locations
                if isfield(data,'spuriousX')
                    X = data.spuriousX;
                    Y = data.spuriousY;
                end
                % get list of CC failure locations
                fail = logical(bitget(data.Status,2));
                X = [X;data.X(fail)]; %#ok<AGROW>
                Y = [Y;data.Y(fail)]; %#ok<AGROW>
                % show results
                plot(X,Y,options.lineSpec,'MarkerSize',options.markerSize);
                
            case {'quiver','quivermean','spquiver'}
                % set default option and parse options
                options.lineSpec = '-k';
                options.qScale = -1.5;
                options.qScalePercentile = 0.8;
                options.clipLo = -Inf;
                options.clipHi = +Inf;
                options.subtractU = 0;
                options.subtractV = 0;
                options.selectLo = -Inf;
                options.selectHi = +Inf;
                options.selectStat = 'all';
                options.selectMult = 0;
                options.vecpos = 0.5;
                options.spMinCC = -Inf;
                [options,ki] = parseOptions(inputs,ki+1,options);
                if strcmpi(options.subtractU,'mean'), options.subtractU = meannan(data.U); end
                if strcmpi(options.subtractV,'mean'), options.subtractV = meannan(data.V); end
                % get data
                switch lower(command)
                    case 'quiver'
                        data = pivSlice;
                        if ~isfield(data,'N')
                            data.N = numel(data.X);
                        end
                        X = reshape(data.X,data.N,1);
                        Y = reshape(data.Y,data.N,1);
                        U = reshape(data.U,data.N,1)-options.subtractU;
                        V = reshape(data.V,data.N,1)-options.subtractV;
                        S = reshape(data.Status,data.N,1);
                        dx = data.iaStepX;
                        dy = data.iaStepY;
                        xmin = min(min(data.X))-dx/2;
                        xmax = max(max(data.X))+dx/2;
                        ymin = min(min(data.Y))-dy/2;
                        ymax = max(max(data.Y))+dy/2;
                        
                    case 'quivermean'
                        data = pivAllSlices;
                        if ~isfield(data,'N')
                            data.N = numel(data.X);
                        end
                        X = reshape(data.X,data.N,1);
                        Y = reshape(data.Y,data.N,1);
                        if isfield(data,'Umean')
                            U = reshape(data.Umean,data.N,1)-options.subtractU;
                            V = reshape(data.Vmean,data.N,1)-options.subtractV;
                        else
                            U = reshape(mean(data.U,3),data.N,1)-options.subtractU;
                            V = reshape(mean(data.V,3),data.N,1)-options.subtractV;
                        end
                        S = zeros(data.N,1);
                        dx = data.iaStepX;
                        dy = data.iaStepY;
                        xmin = min(min(data.X))-dx/2;
                        xmax = max(max(data.X))+dx/2;
                        ymin = min(min(data.Y))-dy/2;
                        ymax = max(max(data.Y))+dy/2;
                    case 'spquiver'
                        data = pivAllSlices;
                        data.N = numel(data.spX);
                        X = reshape(data.spX,data.N,1);
                        Y = reshape(data.spY,data.N,1);
                        Ccmax = reshape(data.spCCmax,data.N,1);
                        U = reshape(mean(data.spUfit,3),data.N,1)-options.subtractU;
                        V = reshape(mean(data.spVfit,3),data.N,1)-options.subtractV;
                        S = zeros(data.N,1);
                        dx = data.spX(2,2)-data.spX(1,1);
                        dy = data.spY(2,2)-data.spY(1,1);
                        xmin = min(min(data.spX))-dx/2;
                        xmax = max(max(data.spX))+dx/2;
                        ymin = min(min(data.spY))-dy/2;
                        ymax = max(max(data.spY))+dy/2;
                        auxNOK = logical(Ccmax<options.spMinCC);
                        U(auxNOK) = NaN;
                        V(auxNOK) = NaN;
                end
                if isfield(data,'multiplicator')
                    M = reshape(data.multiplicator(:,:,1),data.N,1);
                else
                    M = zeros(data.N,1);
                end

                % clip velocity magnitude
                Umag = sqrt(U.^2+V.^2);
                localScale = ones(data.N,1);
                aux = logical(Umag>options.clipHi);
                localScale(aux) = options.clipHi./Umag(aux);
                aux = logical(Umag<options.clipLo);
                localScale(aux) = options.clipLo./Umag(aux);
                % compute the scale
                auxUmag = sort(Umag(~isnan(Umag)));
                auxUmag = auxUmag(round(options.qScalePercentile*numel(auxUmag)));
                if isnan(options.qScale) || options.qScale == 0
                    options.qScale = 2 * sqrt(dx.^2+dy.^2) ./ auxUmag;
                elseif options.qScale < 0
                    options.qScale = -options.qScale * sqrt(dx.^2+dy.^2) ./ auxUmag;
                end
                % select data - Umag
                ok = logical(((~isnan(U)) | (~isnan(V))) & ...
                    (Umag>options.selectLo) & (Umag < options.selectHi));
                % select data - status
                switch lower(options.selectStat)
                    case 'all'
                    case 'valid'
                        ok = ok & ~(logical(bitget(S,5)) | logical(bitget(S,8)));
                    case 'replaced'
                        ok = ok & (logical(bitget(S,5)) | logical(bitget(S,5)));
                    case 'ccfailed'
                        ok = ok & ~(logical(bitget(S,3)) | logical(bitget(S,2)));
                    case 'invalid'
                        ok = ok & (logical(bitget(S,4))) | (logical(bitget(S,7)));
                end
                % select data - multiplicator
                if options.selectMult ~= 0
                    ok = ok & logical(M==options.selectMult);
                end
                % apply scale
                U = U.*localScale*options.qScale;
                V = V.*localScale*options.qScale;
                % change vector position
                X = X - options.vecpos*U;
                Y = Y - options.vecpos*V;
                % apply selection
                X = X(ok);
                Y = Y(ok);
                U = U(ok);
                V = V(ok);
                % select only vectors in frame
                if isfield(pivData,'imSizeX')
                    inFrame = logical(...
                        (X>=xmin ) & (Y>=ymin ) & (X<=xmax) & (Y<=ymax) & ...
                        (X+U>=xmin) & (Y+V>=ymin) & (X+U<=xmax) & (Y+V<=ymax));
                else
                    inFrame = true(size(X));
                end
                % plot quiver
                quiver(X(inFrame), Y(inFrame), U(inFrame), V(inFrame), 0, options.lineSpec);
                axis equal;
                set(gca,'YDir','reverse');
                
            case 'title'
                try
                    titletext = inputs{ki+1};
                    ki = ki+1;
                    title(titletext);
                catch
                    disp('Warning (pivQuiver.m): Title badly specified; ignoring "title" options.');
                end
                
            otherwise
                fprintf('Warning (pivQuiver.m): Unable to parse input "%s". Ignoring it.\n',inputs{ki});
        end
    else
        fprintf('Warning (pivQuiver.m): Unable to parse input %dth input. Ignoring it.\n',ki);
    end
    ki = ki+1;
end

if holdStatus, hold on; else hold off; end
drawnow;
end


%% Local functions

function [options, kout] = parseOptions(input,kin,defaults)
% parseInput - extract commands from input cell
names = fieldnames(defaults);
options = defaults;
success = numel(input)-kin+1>=2;
while success
    try
        command = input{kin};
        value = input{kin+1};
        success = false;
        for jj = 1:numel(names)
            if strcmpi(command,names{jj})
                success = true;
                options.(names{jj}) = value;
                kin = kin + 2;
                break;
            end
        end
    catch
        success = false;
    end
end
kout = kin-1;
end


function [cm] = vorticityColorMap
% color map used for vorticity display
cm = ones(256,3);
cm(129:end,2) = (1:-1/127:0)';
cm(129:end,3) = (1:-1/127:0)';
cm(1:128,1) = (0:1/127:1)';
cm(1:128,2) = (0:1/127:1)';
end


function [out] = meannan(in)
% compute mean of in, ignoring all NaN'a in it
in = reshape(in,numel(in),1);
in = in(~isnan(in));
out = mean(in);
end

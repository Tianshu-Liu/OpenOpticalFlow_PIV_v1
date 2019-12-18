function [pivDataOut] = pivCorrector(pivData,pivData0,pivPar)
% pivCorrector - correct the changes to the velocity field for attenuation

pivDataOut = pivData;

% exit, if corrector is not required, or if initial velocity data are not available
if (~isfield(pivPar,'crAmount')) || ...
        (pivPar.crAmount == 0) || ...
        (~isfield(pivData0,'U')) || ...
        (~isfield(pivData0,'V')) || ...
        (~isfield(pivData0,'X')) || ...
        (~isfield(pivData0,'Y'))
    return;
end

% simplify names
X = pivData.X;
Y = pivData.Y;
U = pivData.U;
V = pivData.V;
W0 = pivData.ccW;

% create interpolators for initial and actual deformation fields
U0estimator = griddedInterpolant(pivData0.X',pivData0.Y',inpaint_nans(pivData0.U'),'spline');
V0estimator = griddedInterpolant(pivData0.X',pivData0.Y',inpaint_nans(pivData0.V'),'spline');
Uestimator = griddedInterpolant(pivData.X',pivData.Y',inpaint_nans(pivData.U'),'spline');
Vestimator = griddedInterpolant(pivData.X',pivData.Y',inpaint_nans(pivData.V'),'spline');

% compute for all image pixels
[XX,YY] = ndgrid(1:pivData.imSizeX,1:pivData.imSizeY);
U0im = U0estimator(XX,YY)';
V0im = V0estimator(XX,YY)';
Uim = Uestimator(XX,YY)';
Vim = Vestimator(XX,YY)';

% get image mask
auxM1 = pivData.imMaskArray1;
auxM2 = pivData.imMaskArray2;
if numel(auxM1)==0, auxM1 = ones(size(XX)); end
if numel(auxM2)==0, auxM2 = ones(size(XX)); end
M = auxM1.*auxM2;

    % interpolate initial guess for final grid
U0 = U0estimator(X',Y')';  
V0 = V0estimator(X',Y')';

% loop over all IAs and correct the velocity
for kx = 1:size(X,2)
    for ky = 1:size(Y,1)
        % get mask of the current IA
        auxM = M(...
            Y(ky,kx)-(pivPar.iaSizeY-1)/2+1:Y(ky,kx)+(pivPar.iaSizeY-1)/2+1, ...
            X(ky,kx)-(pivPar.iaSizeX-1)/2+1:X(ky,kx)+(pivPar.iaSizeX-1)/2+1 ...
            );
        % get actual weigting function, respecting the mask
        W = W0.*auxM;
        % apply corrector only if less than 75% of "weighted" pixels are masked
        if sum(sum(W)) > 0.25*sum(sum(W0))
            % normalize weigting function
            W = W/sum(sum(W));
            % get the weighted averages of velocity estimates in the IA
            U0avg = sum(sum(U0im(...
                Y(ky,kx)-(pivPar.iaSizeY-1)/2+1:Y(ky,kx)+(pivPar.iaSizeY-1)/2+1, ...
                X(ky,kx)-(pivPar.iaSizeX-1)/2+1:X(ky,kx)+(pivPar.iaSizeX-1)/2+1 ...
                ).*W));
            V0avg = sum(sum(V0im(...
                Y(ky,kx)-(pivPar.iaSizeY-1)/2+1:Y(ky,kx)+(pivPar.iaSizeY-1)/2+1, ...
                X(ky,kx)-(pivPar.iaSizeX-1)/2+1:X(ky,kx)+(pivPar.iaSizeX-1)/2+1 ...
                ).*W));
            Uavg = sum(sum(Uim(...
                Y(ky,kx)-(pivPar.iaSizeY-1)/2+1:Y(ky,kx)+(pivPar.iaSizeY-1)/2+1, ...
                X(ky,kx)-(pivPar.iaSizeX-1)/2+1:X(ky,kx)+(pivPar.iaSizeX-1)/2+1 ...
                ).*W));
            Vavg = sum(sum(Vim(...
                Y(ky,kx)-(pivPar.iaSizeY-1)/2+1:Y(ky,kx)+(pivPar.iaSizeY-1)/2+1, ...
                X(ky,kx)-(pivPar.iaSizeX-1)/2+1:X(ky,kx)+(pivPar.iaSizeX-1)/2+1 ...
                ).*W));
            % correct the velocity change
            pivDataOut.U(ky,kx) = U(ky,kx) + pivPar.crAmount*((U(ky,kx)-U0(ky,kx))-(Uavg-U0avg));
            pivDataOut.V(ky,kx) = V(ky,kx) + pivPar.crAmount*((V(ky,kx)-V0(ky,kx))-(Vavg-V0avg));
        end            
    end
end
end
%% kernelSmoother
% Implements kernel smoothing for non-parametric regression.
%
%   Input arguments:
%       X:  X-Samples
%       Y:  g(X)-Samples
%       options: MCAssetPricing options object.
%       X0: Evaluation points of regressed function (default: X)
%   Output arguments:
%       YEstimates: Estimates of f(X0) based on kernel smoothing. 
function YEstimates = kernelSmoother(X, Y, options, X0) 
    if nargin < 4
        X0 = X;
    end
    
    % Always reshape input arrays to matrices.
    S = size(X0);
    M = S(1);
    D = prod(S(2:end));
    X0 = reshape(X0,[M,D]);
    X = reshape(X,[size(X,1),D]);
    YEstimates = zeros(M,1);
    
    % Kernel Selection. Default: Gaussian kernel.
    K = options.smoothingKernel;

    % For each estimation point.
    for m = 1:M
        % Evaluate kernel function for each X(i) and compute weight.
        W = K(X0(m,:),X,options);

        % Fitting constant function (Kernel Smoother).
        if (options.basisDegree == 0)
            YEstimates(m) = (W'*Y)/sum(W);
            
        % Local regression onto basis.
        else
            % Compute weighted model matrix. Default: Polynomial
            A = diag(sqrt(W)) * options.computeBasisVector(X,options);

            % Fit to model.
            c = A\(diag(sqrt(W))*Y);

            % Estimate Y(X0(m,:)).
            YEstimates(m) = options.computeBasisVector(X0(m,:),options) * c;
        end
    end
end

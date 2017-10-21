%% generateJumpDiffusionPaths
% Generates sample paths of a D-dimensional geometric brownian motion with jumps.
%
%   Input arguments:
%       options: MCAssetPricing options object.
%		volatility: [1,M] matrix to be used in place of sigma (for stochastic volatility)
%   Output arguments:
%       S: [N,M,D] matrix representing N samples of M timesteps of a D
%           dimensional geometric brownian motion.
function S = generateJumpDiffusionPaths(options, volatility)

    % Dimensions.
    [N, M] = deal(options.numberOfSamples, options.numberOfTimesteps);

    % Black-Scholes parameters.
    [r, dividends, sigma, S0, Corr] = deal(options.r, options.dividends, options.sigma, options.S0, options.correlationMatrix);
    
    if (nargin == 1) 
        volatility = repmat(sigma, [1 M]);
    end
    
    % Jump process parameters.
    [phi, lambda] = deal(options.phi, options.lambda);
    
    % Generator parameters.
    dt = options.dt;

    % Initialize matrix.
    D = numel(S0);
    S = zeros(N,M,D);
    S(:,1,:) = repmat(S0,[N,1]);

    % Cholesky decomposition.
    Rho = chol(Corr);
    
    % Propagate.
    for t = 2:M
        sigma = volatility(t);
        
        % Uncorrelated Normal samples. Default: randn
        B = options.generateNormalSamples(N,D);

        % Correlate Normal samples.
        W = B * Rho;
        
        % Scale Normal samples.
        W = W * diag(sigma);

        % Propagation step.
        S(:,t,:) = squeeze(S(:,t-1,:)) .* exp(repmat(r(t-1)-dividends-.5*(sigma.^2), [N 1])*dt + sqrt(dt).*W);
        
        % Jump step.
        if any(lambda ~= 0 && phi ~= 0)
            S(:,t,:) = squeeze(S(:,t,:)) .* exp(log(1 + phi).*poissrnd(repmat(lambda.*dt,[N,1])));
        end
    end
end
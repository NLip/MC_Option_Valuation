%% generateHestonPaths
% Generates sample paths of a 1-dimensional Heston path (stochastic volatility).
%   Input arguments:
%       options: MCAssetPricing options object.
%   Output arguments:
%       S: [N,M] matrix representing N samples of M timesteps of a 1
%           dimensional geometric brownian motion.
function S = generateHestonPaths(options)
   	% Dimensions.
    [N, M] = deal(options.numberOfSamples, options.numberOfTimesteps);

    % Parameters.
    [kappa,theta, eta, sigma, T] = deal(options.kappa, options.theta, options.eta, options.sigma, options.T);
    
    a = @(t,v) kappa*(theta - v);
    b = @(t,v) (eta * sqrt(v));
    dt = (T/M)/240; % Simulate in 6 minute intervals
    volatility = sqrt(eulerMethod(a,b,sigma,M,T,N,dt));
    
    S = generateJumpDiffusionPaths(options, volatility);
    
end
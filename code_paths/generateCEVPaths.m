%% generateCEVPaths
% Generates sample paths of a 1-dimensional of the solution to dSt = St (rdt + sigma S^(1-gamma)dWt).
%
%   Input arguments:
%       options: MCAssetPricing options object.
%   Output arguments:
%       S: [N,M] matrix representing N samples of M timesteps of a 1
%           dimensional CEV path.
function S = generateCEVPaths(options)
   	% Dimensions.
    [N, M] = deal(options.numberOfSamples, options.numberOfTimesteps);

    % Parameters.
    [r, dividends, sigma, S0, T, gamma] = deal(options.r, options.dividends, options.sigma, options.S0, options.T, options.gamma);
    
    a = @(t,s) (r(1) - dividends)*s; % TODO: Allow changing interest (not have r(1))
    b = @(t,s) sigma*s.^gamma;
    dt = (T/M)/240; % Simulate in 6 minute intervals
    S = eulerMethod(a,b,S0,M,T,N,dt);
    
end
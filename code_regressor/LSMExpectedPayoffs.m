%% LSMExpectedPayoffs
% Implements least squares monte carlo computation of expected payoffs
%
%   Input arguments:
%       St: Sample path at a given timestep.
%       discountedCashFlows: TODO
%       options: MCAssetPricing options object.
%   Output arguments:
%       
function expectedPayoffs = LSMExpectedPayoffs(St, discountedCashFlows, options)
    
    % Computes model matrix. Default: Polynomial
    A = options.computeBasisVector(St,options);

    % Least squares fit.
    coeffs = A\discountedCashFlows;
    
    % Apply model.
    expectedPayoffs = A*coeffs;
end
%% weightedLaguerrePolynomial
% Evaluates weighted Laguerre basis in points St (in one dimension)
%
%   Input arguments:
%		X: Evaluation points
%       options: MCAssetPricing options object.
%   Output arguments:
%       A: [N,basisDegree] matrix representing evaluation of basis function for each of the N samples
function A = weightedLaguerrePolynomial(XX,options)
    D = options.basisDegree;
    X = XX./options.S0;
    
    A = zeros(numel(X),D + 1);
    A(:,1) = 1;
    A(:,2) = 1 - X;
    for k = 2:D
        A(:,1+k) = ((2*k - 1 - X).*A(:,k) - (k-1).*A(:,k-1))/(k);
    end
    A = repmat(exp(-X./2),[1 D+1]).*A;
end
        
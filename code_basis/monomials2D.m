%% monomials2D
% Evaluates 2D monomial basis in points St
%
%   Input arguments:
%		St: Evaluation points
%       options: MCAssetPricing options object.
%   Output arguments:
%       S: [N,nchoosek(basisDegree + D,D)] matrix representing evaluation of basis function for each of the N samples
function result = monomials2D(St, options)
    [N,D] = size(St);
    assert(D == 2);
    x = St(:,1);
    y = St(:,2);
    
    basisDegree = options.basisDegree;
    offset = 0;
    result = zeros(N,nchoosek(basisDegree + D,D));
    for order = 0:basisDegree
        numberOfMonomials = order + 1;
        xx = fliplr(cumprod([ones(size(x)) repmat(x,[1 order])],2));
        yy = cumprod([ones(size(y)) repmat(y,[1 order])],2);
        result(:,offset + (1:numberOfMonomials)) = xx.*yy;
        offset = offset + numberOfMonomials;
    end
end
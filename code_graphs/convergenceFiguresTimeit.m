function data = convergenceFiguresTimeit(data)
    
    options = MCAssetPricingOptions('seed',0.42);
    if ~isfield(data,'reference')
        data.reference = AmericanOptFD(options.S0,options.K,options.r(1),options.T,options.sigma,2000,2000);
    end
    
    P = 20; Nmin = 3; Nmax = 150000;
    N = floor(linspace(Nmin,Nmax,P));


    makePlot(...
        'laguerre3',...
        MCAssetPricingOptions(options,'computeBasisVector',@weightedLaguerrePolynomial,'basisDegree',3),...
        N, ...
        data.reference);
%     
%         makePlot(...
%         'laguerre2',...
%         MCAssetPricingOptions(options,'computeBasisVector',@weightedLaguerrePolynomial,'basisDegree',2),...
%         N, ...
%         data.reference);
%     
%         makePlot(...
%         'laguerre4',...
%         MCAssetPricingOptions(options,'computeBasisVector',@weightedLaguerrePolynomial,'basisDegree',4),...
%         N, ...
%         data.reference);
%     
%         makePlot(...
%         'laguerre6',...
%         MCAssetPricingOptions(options,'computeBasisVector',@weightedLaguerrePolynomial,'basisDegree',6),...
%         N, ...
%         data.reference);
%     
%         makePlot(...
%         'laguerre12',...
%         MCAssetPricingOptions(options,'computeBasisVector',@weightedLaguerrePolynomial,'basisDegree',12),...
%         N, ...
%         data.reference);
% 
%         makePlot(...
%         'heston',...
%         MCAssetPricingOptions(options,'generatePaths',@generateHestonPaths,'sigma',options.sigma^2, 'kappa',0.2,'eta',0.4),...
%         N, ...
%         data.reference);
% 
%         makePlot(...
%         'CEV',...
%         MCAssetPricingOptions(options,'generatePaths',@generateCEVPaths,'sigma',options.S0^(1-options.gamma)*options.sigma),...
%         N, ...
%         data.reference);
    
%         makePlot(...
%         'polynomial1',...
%         MCAssetPricingOptions(options,'basisDegree',1),...
%         N, ...
%         data.reference)
    
%             makePlot(...
%         'polynomial2',...
%         MCAssetPricingOptions(options,'basisDegree',2),...
%         N, ...
%         data.reference)
%     
%             makePlot(...
%         'polynomial3',...
%         MCAssetPricingOptions(options,'basisDegree',3),...
%         N, ...
%         data.reference)
%     
%             makePlot(...
%         'polynomial4',...
%         MCAssetPricingOptions(options,'basisDegree',4),...
%         N, ...
%         data.reference)
%     
%             makePlot(...
%         'polynomial6',...
%         MCAssetPricingOptions(options,'basisDegree',6),...
%         N, ...
%         data.reference)

%         doScrambling = false;
%         getSobolPoints = initSobol(doScrambling);    
%         makePlot(...
%         'sobol',...
%         MCAssetPricingOptions(options,'generateNormalSamples', @(N,M) norminv(getSobolPoints(N))),...
%         N, ...
%         data.reference)
%     
%         doScrambling = true;
%         getSobolPoints = initSobol(doScrambling);    
%         makePlot(...
%         'sobol',...
%         MCAssetPricingOptions(options,'generateNormalSamples', @(N,M) norminv(getSobolPoints(N))),...
%         N, ...
%         data.reference)
% 
%         D = 2;
%         makePlot(...
%         '2dimensional',...
%         MCAssetPricingOptions(options,'computeBasisVector',@monomials2D,'basisDegree',2,...
%             'S0',repmat(options.S0,[1,D]), ...
%             'dividends',repmat(options.dividends,[1,D]),...
%             'sigma',repmat(options.sigma,[1,D])*sqrt(D),...
%             'correlationMatrix',eye(D), ...
%             'computePayoffs',@(S,t,options) (max(options.K-mean(S(:,t,:),3),0))...
%         ),...
%         N, data.reference)
%     
% 
%     
%         tau = 5;
%         makePlot(...
%         'timeDependency',...
        MCAssetPricingOptions(options,'basisDegree',3,'isExercisable',@(t,options) (t > tau),'getMemory',@(t) (t-tau):t,  'computePayoffs',@(S,t,options) (max(options.K - sum(S(:,((t-tau):t)),2)/(tau + 1),0)) 
%         N)
% 
%         makePlot(...
%         'marketsCrashAndBurn',...
%         MCAssetPricingOptions(options,...
%             'phi',-0.6, ...
%             'lambda',1 ...
%         ),...
%         N)
%     
%     
%     
%     makePlot(...
%         'kernelSmoother',...
%         MCAssetPricingOptions(options,'computeExpectedPayoffs', @kernelSmoother,'basisDegree',0),...
%         N, ...
%         data.reference);
%     
%         D = 3;
%         makePlot(...
%         '3dimensional',...
%         MCAssetPricingOptions(options,'computeExpectedPayoff',@kernelSmoother,'basisDegree',0,...
%             'S0',repmat(options.S0,[1,D]), ...
%             'dividends',repmat(options.dividends,[1,D]),...
%             'sigma',repmat(options.sigma,[1,D])*sqrt(D),...
%             'correlationMatrix',eye(D), ...
%             'computePayoffs',@(S,t,options) (max(options.K-mean(S(:,t,:),3),0))...
%         ),...
%         N, data.reference)
end

function makePlot(name, options, N, referenceSolution)
    figure()
    value = zeros(size(N));
    try
        for i = 1:numel(N);
            [value(i), ~, ~, exerciseTimes] = MCAssetPricing(options,'numberOfSamples',N(i));
            disp(i);
        end
        if (nargin > 3) 
            subplot(1,3,1)
            plot([N(1),N(end)],repmat(referenceSolution,[1,2]),'g-')
            hold on;
        else
            subplot(1,2,1)
        end

        plot(N,value,'b*');
        if (nargin > 3)
            subplot(1,3,2)
        else
            subplot(1,2,2)
        end
        histogram(exerciseTimes,1:options.numberOfTimesteps);
        
        if (nargin > 3)
            subplot(1,3,3)
            error = abs(value-referenceSolution);
            semilogy(N,error,'r*');
        end
        savefig(name);
        matlab2tikz(strcat(name,'.tex'));
    catch error
        disp(name)
        disp(error)
    end
end

function f = initSobol(doScramble)
    f = @getSobolPoints;
    P = 0;
    
    function X = getSobolPoints(N) 
        S = sobolset(1,'Skip',1e3 + P);
        if doScramble
            S = scramble(S,'MatousekAffineOwen');
        end
        P = P + N;
        X = net(S,N);
    end
end
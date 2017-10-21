function tests(testnumbers)
    M = 2; % Total number of tests
    if nargin == 0
        testnumbers = 1:M;
    end
    
    for test = testnumbers
        fprintf('Running test number: %i\n',test)
        switch (test)
            case 1
                figure()
                fprintf('\tAmerican Put option, default settings:')
                [value,~,~,~,S,~] = MCAssetPricing();
                value
                semilogy(S')
                
            case 2
                % American Call option.
                fprintf('\tAmerican Call option price:')
                MCAssetPricing(...
                    'computePayoffs', @(S,t,options) (max(S(:,t) - options.K,0)), ...
                    'r',0.04,'K', 80, 'S0', 100,'seed',42,'basisDegree',3)

            case 4
                fprintf('\t"Reference Solution": %g\n', ...
                     MCAssetPricing('numberOfSamples',1e6,'basisDegree',3 ...
                ));
            case 5
                options = MCAssetPricingOptions('numberOfSamples',10000,'numberOfTimesteps',3, ...
                    'correlationMatrix',[1 0.5; 0.5 1],'sigma',[0.3 0.4], ...
                    'dividends',[0 0],'S0',[80, 80],'seed',0,...
                    'computeBasisVector',@monomials2D, ...
                    'computePayoffs', @(S,t,options) max(options.K  - min(S(:,t,:),[],3),0)...
                );
                MCAssetPricing(options)
            case 6
                MCAssetPricing('numberOfSamples',10000,'computeExpectedPayoffs', @kernelSmoother, 'basisDegree',0)
            case 7
                % Multidimensional option with memory.
                tau = 2;
                options = struct( ...
                    'seed',                     42, ...
                    'basisDegree',              0, ...
                    'generatePaths',            @generateBlackScholesPaths, ...
                    'S0',                       [80,80], ...
                    'K',                        100, ...
                    'r',                        0.04, ...
                    'dividends',                [0 0], ...
                    'sigma',                    [0.2,0.2], ...
                    'correlationMatrix',        eye(2), ...
                    'computePayoffs',           @(S,t,options) (max(options.K - sum(S(:,((t-tau):t),1),2)/(tau + 1),0)), ...
                    'computeExpectedPayoffs',   @kernelSmoother, ...
                    'isExercisable',            @(t,options) (t > tau), ...
                    'getMemory',                @(t) (t-tau):t ...
                );
                MCAssetPricing(options)
            case 8
                M = 100;
                val = zeros(1,M);
                for i = 1:M
                   val(i) = MCAssetPricing('numberOfSamples',1000);
                end
                plot(val)
            case 9
                S0 = 80;
                K = 100;
                r = 0.04;
                sigma = 0.2;
                T = 1;
                D = 20;
                value = zeros(1,D);
                Nmax = 10000;
                N = (10:floor(Nmax/D):Nmax);
                M = N;
                for i = 1:D
                    value(i) = AmericanOptFD(S0,K,r,T,sigma,N(i),M(i));
                end
                semilogx(value);
            case 10
                
                S0 = 80;
                K = 100;
                r = 0.04;
                sigma = 0.0;
                T = 1;
                M = 10000;
                AmericanOptFD(S0,K,r,T,sigma,M,M)
                MCAssetPricing('S0',S0,'K',K,'r',r,'sigma',sigma,'T',T,'numberOfSamples',M)
            case 11
                MCAssetPricing('numberOfSamples',1,'S0',1,'sigma',0,'r',0,'lambda',4)
                
        end
    end
end

function testMCAssetPricing(testnumbers)
    M = 12; % Total number of tests
    if nargin == 0
        testnumbers = 1:M;
    end
    
    for test = testnumbers
        fprintf('Running test number: %i\n',test)
        switch (test)
            case 1
                figure()
                fprintf('\tAmerican Put option, default settings:')
                [val,~,~,~,S,opt] = MCAssetPricing('seed',0);
                val
                semilogy(S')
                
            case 2
                % American Call option.
                fprintf('\tAmerican Call option price:')
                MCAssetPricing(...
                    'computePayoffs', @(S,t,options) (max(S(:,t) - options.K,0)), ...
                    'r',0.04,'K', 80, 'S0', 100,'seed',42,'basisDegree',3)
            case 3
                figure()
                fprintf('Sobol sequence')

                doScrambling = false;
                getSobolPoints = initSobol(doScrambling);
                [val,~,~,~,S] = MCAssetPricing(...
                    'generateNormalSamples', @(N,M) norminv(getSobolPoints(N)), ...
                    'seed', 0.42, 'numberOfSamples',9871 ...
                );
                val
                subplot(2,2,3)
                semilogy(S')
                
                getSobolPoints = initSobol(doScrambling);
                [val,~,~,~,S] = MCAssetPricing(...
                    'generateNormalSamples', @(N,M) norminv(getSobolPoints(N)), ...
                    'seed', 0.42, 'numberOfSamples', 9878 ...
                );
                val
                subplot(2,2,4)
                semilogy(S')
                
                getSobolPoints = initSobol(doScrambling);
                [val,~,~,~,S,options] = MCAssetPricing(...
                    'generateNormalSamples', @(N,M) norminv(getSobolPoints(N)), ...
                    'seed', 0.42, 'numberOfSamples',1000 ...
                );
                val
                subplot(2,2,1)
                semilogy(S')
                
                getSobolPoints = initSobol(doScrambling);
                [val,~,~,~,S,options] = MCAssetPricing(...
                    'generateNormalSamples', @(N,M) norminv(getSobolPoints(N)), ...
                    'seed', 0.42, 'numberOfSamples',987 ...
                );
                val
                subplot(2,2,2)
                semilogy(S')
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
            case 12
                doScrambling = false;
                getSobolPoints = initSobol(doScrambling);
                [val,~,~,~,SY,options] = MCAssetPricing(...
                    'generateNormalSamples', @(N,M) norminv(getSobolPoints(N)), ...
                    'seed', 0.42, 'numberOfSamples',987 ...
                );
            
                N = options.numberOfSamples;
                X = randn(N,1);
                Y = options.generateNormalSamples(N,1);
                
                doScrambling = true;
                getSobolPoints = initSobol(doScrambling);
                [val,~,~,~,SZ,options] = MCAssetPricing(...
                    'generateNormalSamples', @(N,M) norminv(getSobolPoints(N)), ...
                    'seed', 0.42, 'numberOfSamples',987 ...
                );
            
                Z = options.generateNormalSamples(N,1);
                
                subplot(3,3,1)
                qqplot(X)
                subplot(3,3,2)
                qqplot(Y)
                subplot(3,3,3)
                qqplot(Z)
                
                subplot(3,3,4)
                hist(X)
                
                subplot(3,3,5)
                hist(Y)
                
                subplot(3,3,6)
                hist(Z)

                subplot(3,3,7)
                [val,~,~,~,SX,opt] = MCAssetPricing('seed',0.42, 'numberOfSamples',987);
                semilogy(SX')
                
                subplot(3,3,8)
                semilogy(SY')
                
                subplot(3,3,9)
                semilogy(SZ')
                
                varX = var(X)
                varY = var(Y)
                varY = var(Z)
                diffMeanY = mean(X) - mean(Y)
                diffMeanZ = mean(X) - mean(Z)
                
            case 13
                figure()
                [val,~,~,~,S,opt] = MCAssetPricing('lambda',0.2, 'phi',-0.6);
                val
                semilogy(S')
        end
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

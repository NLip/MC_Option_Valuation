function options = MCAssetPricingOptions(varargin)

    % -------------------------- Default Options ----------------------------- %
    
    rng shuffle % To ensure that defaults.seed is chosen differently everytime. 
    defaults = struct( ...
        'seed',                     randi(intmax), ... % with sobol try different seeds
        'numberOfSamples',          1000, ... % TODO: Convergence
        'computeBasisVector',       @(St, options) cumprod([ones(size(St)) repmat(St,[1 options.basisDegree])],2), ... % also try laguerre
        'basisDegree',              3, ... % Maybe try [1, 3, 6, 12]
        'dt',                       1/50, ...
        'T',                        1, ...
        'numberOfTimesteps',        50, ...
        'generatePaths',            @generateJumpDiffusionPaths, ...
        'S0',                       80, ...
        'K',                        100, ...
        'r',                        0.01, ...
        'dividends',                0, ...
        'sigma',                    0.25, ...
        'correlationMatrix',        1, ...
        'generateNormalSamples',    @randn, ... % Try out sobol sequences
        'phi',                      -0.1, ...
        'lambda',                   0, ...
        'computePayoffs',           @(S,t,options) (max(options.K-S(:,t,1),0)), ... % TODO: Other options (asian)
        'isInTheMoney',             @(S,t,options) (options.computePayoffs(S,t,options) > 0), ...
        'isExercisable',            @(t,options) true, ... % TODO: other options (bermuda)
        'computeExpectedPayoffs',   @LSMExpectedPayoffs, ...  % maybe change to kernel smoother (set basisDegree to 0)
        'computeWhoExercises',      @whoExercisesExpectedPayoff, ... 
        'getMemory',                @(t) t, ...
        'smoothingBandwidth',       3, ...
        'smoothingKernel',          @(X0,X,options) exp(-sum((repmat(X0,[size(X,1),1]) - X).^2./(2*options.smoothingBandwidth^2),2)), ...
        'kappa', 					2, ...
        'theta', 					0.25^2, ...
        'eta', 						0.02, ...
        'gamma',                    1.2 ...
    );

    % --------------------------- Setup Parser ------------------------------- %
    
    % Custom conditions.
    isfunction = @(x) isa(x,'function_handle');
    ispositiveInteger = @(x) isnumeric(x) && (x > 0) && (x == floor(x));
    isnonnegativeInteger = @(x) isnumeric(x) && (x >= 0) && (x == floor(x));
    ispositiveNumber = @(x) isnumeric(x) && (x > 0);
    iscorrelationMatrix = @(x) isnumeric(x) && all(diag(x) == 1) && ispositiveDefinite(x) && issymmetric(x); 
    
    % Initialize parser.
    p = inputParser();

    % Product Parameters.
    addParameter(p,'computePayoffs',defaults.computePayoffs,isfunction);
    addParameter(p,'isInTheMoney',defaults.isInTheMoney,isfunction);
    addParameter(p,'isExercisable',defaults.isExercisable,@(x) isfunction(x) || isnumeric(x));
    addParameter(p,'K',defaults.K,@isnumeric);
    addParameter(p,'dividends',defaults.dividends,@isnumeric);
    addParameter(p,'getMemory',defaults.getMemory,isfunction);

    % Simulation Parameters.
    addParameter(p,'seed',defaults.seed);        	
    addParameter(p,'generateNormalSamples',defaults.generateNormalSamples,isfunction);
    addParameter(p,'computeBasisVector',defaults.computeBasisVector,isfunction);
    addParameter(p,'computeExpectedPayoffs',defaults.computeExpectedPayoffs,isfunction);
    addParameter(p,'computeWhoExercises',defaults.computeWhoExercises,isfunction);
    addParameter(p,'smoothingBandwidth',defaults.smoothingBandwidth,@isnumeric);
    addParameter(p,'smoothingKernel',defaults.smoothingKernel,isfunction);
    addParameter(p,'basisDegree',defaults.basisDegree, isnonnegativeInteger);
    addParameter(p,'numberOfSamples',defaults.numberOfSamples, ispositiveInteger);
    addParameter(p,'numberOfTimesteps',defaults.numberOfTimesteps, ispositiveInteger);

    % Path Generation Parameters.
    addParameter(p,'T',defaults.T,ispositiveNumber);
    addParameter(p,'dt',defaults.dt,ispositiveNumber);
    addParameter(p,'S0',defaults.S0,@isnumeric);
    addParameter(p,'r',defaults.r,@isnumeric);
    addParameter(p,'sigma',defaults.sigma,@isnumeric);
    addParameter(p,'phi',defaults.phi,@isnumeric);
    addParameter(p,'lambda',defaults.lambda,@isnumeric);
    addParameter(p,'kappa',defaults.kappa,@isnumeric);
    addParameter(p,'theta',defaults.theta,@isnumeric);
    addParameter(p,'eta',defaults.eta,@isnumeric);
    addParameter(p,'gamma',defaults.gamma,@isnumeric);
    addParameter(p,'generatePaths',defaults.generatePaths,isfunction);
    addParameter(p,'correlationMatrix',defaults.correlationMatrix,iscorrelationMatrix);

    % Parse input.
    parse(p,varargin{:});
    options = p.Results;
    
    % ----------------------- Check Dimensions ------------------------------- %
    [S0, sigma, correlationMatrix, dividends, lambda, phi] = ...
        deal(options.S0, options.sigma, options.correlationMatrix, options.dividends, options.lambda, options.phi);
    D = numel(S0);
    assert(all(size(S0) == [1,D]), 'S0 must have size [1,D]');
    assert(all(size(sigma) == [1,D]), 'sigma must have size [1,D]');
    assert(all(size(dividends) == [1,D]), 'dividends must have size [1,D]');
    assert(all(size(correlationMatrix) == [D,D]), 'correlationMatrix must have size [D,D]');
    if any(lambda ~= 0 && phi ~= 0)
        assert(all(size(lambda) == [1,D]), 'lambda must have size [1,D]');
        assert(all(size(phi) == [1,D]), 'phi must have size [1,D]');
    end
    % ----------------------- Make Consistent -------------------------------- %
    
    if numel(options.r) ~= (options.numberOfTimesteps - 1)
        options.r = repmat(options.r(1),[1 options.numberOfTimesteps-1]);
    end
    
    % Extract relevant options to avoid retyping 'options.'
    [T, dt, N] = deal(options.T, options.dt, options.numberOfTimesteps);
    [Td, dtd] = deal(defaults.T,defaults.dt);

    % If parameters are not consistent
    if (N * dt ~= T)
        % If T was set differently from the default.
        if (T ~= Td)
            % If dt was set.
            if (dt ~= dtd)
                N = T/dt;
            else
                dt = T/N;
            end
        % If dt was set.
        elseif (dt ~= dtd)
            T = N*dt;
        else
            dt = T/N;
        end
    end

    [options.T, options.dt, options.numberOfTimesteps] = deal(T,dt,N);    
end

function flag = ispositiveDefinite(x)
    [~,p] = chol(x);
    flag = p == 0;
end


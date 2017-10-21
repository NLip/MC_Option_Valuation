function [value, stderr, discountedCashFlows, exerciseTimes, S, options] = MCAssetPricing(varargin)

	% --------------------------- Add  Library ------------------------------- %
    addpath('../code_basis','../code_paths','../code_regressor','../code_strategy');
    
    % --------------------------- Handle Input ------------------------------- %

    % Parse input arguments.
    options = MCAssetPricingOptions(varargin{:});
        
    % Set seed.
    rng(options.seed);
    
    % Simulation parameters.
    [dt, r, numberOfTimesteps, numberOfSamples] = deal(options.dt, options.r, options.numberOfTimesteps, options.numberOfSamples);

    % ---------------------------- Simulation -------------------------------- %

    % Generate sample paths. Default: generateBlackScholesPaths
    S = options.generatePaths(options);

    % Compute payoff at maturity. Default: @(S,t,options) (max(options.K-S(:,t),0))
    discountedCashFlows = zeros(numberOfSamples,numberOfTimesteps);
    discountedCashFlows(:,end) = options.computePayoffs(S,numberOfTimesteps,options);
    
    % Prepare container for exercise times.
    exerciseTimes = numberOfTimesteps * ones(numberOfSamples,1);
    
    % Work backwards from maturity.
    for t=(numberOfTimesteps-1):-1:1

        % Discount cash flows.
        discountedCashFlows(:,t) = exp(-r(t)*dt) * discountedCashFlows(:,t+1);

        % Can the option be exercised in this timestep? Default: @(t,options) true
        if (options.isExercisable(t,options))
            
            % Compute payoffs for exercising in this time step. Default: @(S,t,options) (max(options.K-S(:,t),0))
            payoffs = options.computePayoffs(S,t,options);

            % Find samples where exercising is worthwhile.
            whoExercises = options.computeWhoExercises(S, t, discountedCashFlows(:,t), payoffs, options);
            
            % Exercise if worthwhile.
            exerciseTimes(whoExercises) = t;
            discountedCashFlows(whoExercises,t) = payoffs(whoExercises);
        end       
    end

    % ---------------------------- Statistics -------------------------------- %

    value = mean(discountedCashFlows(:,1));
    stderr = std(discountedCashFlows(:,1)   )/sqrt(numberOfSamples);
end




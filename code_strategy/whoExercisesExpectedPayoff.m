%% whoExercisesExpectedPayoff
% Decides what path would exercise based on the estimation of the future expected payoff.
function [whoExercises, whoIsInTheMoney, payoffs, expectedPayoffs] = whoExercisesExpectedPayoff(S, t, discountedCashFlows, payoffs, options)
    
    % Find paths that are in the money. Default: @(S,t,options) (options.computePayoffs(S,t,options) > 0)
    whoIsInTheMoney = find(options.isInTheMoney(S,t,options));

    % If nobody is in the money, nobody exercises.
    if (numel(whoIsInTheMoney) == 0)
       whoExercises = [];
       
    % Estimate expected payoffs conditional on S.
    else
        % Default @(t) t.
        tau = options.getMemory(t);
        
        % Only consider in-the-money paths of the current timestep since
        % expected payoffs are irrelevant for the other paths which will
        % certainly not exercise.
        St = S(whoIsInTheMoney,tau,:);
        payoffs = payoffs(whoIsInTheMoney);
        discountedCashFlows = discountedCashFlows(whoIsInTheMoney);
        
        % Are all samples equal?
        if (all(St == repmat(St(1,:,:),[size(St,1),1,1])))
            % One cannot fit basis with only one sample for explanatory
            % variable. Furthermore, conditional expectetation estimate 
            % degenerates to mean of dependent variable samples in this case.
            expectedPayoffs = mean(discountedCashFlows);
        else
            % Use regression to estimate conditional expectation. Default: LSMExpectedPayoffs
            expectedPayoffs = options.computeExpectedPayoffs(St,discountedCashFlows, options);
        end
        
        % Find paths for which exercising seems better than continuing.
        whoExercises = whoIsInTheMoney(payoffs > expectedPayoffs);       
    end
end

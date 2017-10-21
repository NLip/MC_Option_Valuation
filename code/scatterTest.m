function scatterTest()
    [value, stderr, discountedCashFlows, exerciseTimes, S, options] = MCAssetPricing('seed',0.42,'numberOfSamples',3000,'numberOfTimesteps',100);
    value
    optionsLinear = MCAssetPricingOptions(options,'basisDegree',1);
    optionsKS = MCAssetPricingOptions(options,'basisDegree',0,'computeExpectedPayoff',@kernelSmoother,'smoothingBandwidth',1.5);
    optionsLR = MCAssetPricingOptions(optionsKS, 'basisDegree',1);

    t = floor(options.numberOfTimesteps/2);
    subplot(1,2,1);
    scatter(S(:,t),discountedCashFlows(:,t),'b');
    payoffs = options.computePayoffs(S,t,options);

    [whoExercises, whoIsInTheMoney, ~, expectedPayoffs] = whoExercisesExpectedPayoff(S, t, discountedCashFlows(:,t), payoffs, options);
    subplot(1,2,2);
    hold on;
    scatter(S(whoIsInTheMoney,t),discountedCashFlows(whoIsInTheMoney,t),'b')
    scatter(S(whoIsInTheMoney,t),expectedPayoffs,'r')
    [whoExercises, whoIsInTheMoney, ~, expectedPayoffs] = whoExercisesExpectedPayoff(S, t, discountedCashFlows(:,t), payoffs, optionsLinear);
    scatter(S(whoIsInTheMoney,t),expectedPayoffs,'g')
    [whoExercises, whoIsInTheMoney, ~, expectedPayoffs] = whoExercisesExpectedPayoff(S, t, discountedCashFlows(:,t), payoffs, optionsKS);
    scatter(S(whoIsInTheMoney,t),expectedPayoffs,'m')
    %[whoExercises, whoIsInTheMoney, ~, expectedPayoffs] = whoExercisesExpectedPayoff(S, t, discountedCashFlows, payoffs, optionsLR);
    %scatter(S(whoIsInTheMoney,t),expectedPayoffs,'k')
    
%     MCAssetPricing(optionsLinear)
%     MCAssetPricing(optionsKS)
    

    
    
end

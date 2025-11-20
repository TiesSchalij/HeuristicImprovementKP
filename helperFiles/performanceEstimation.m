function [valueFunctionEstimations, valueFunctionTildeEstimations] = performanceEstimation(policiesToRun, generator, nSamples, nItems, exact, sampleArray)
% performanceEstimation estimates the performance of different policies on the same sample
% Input:
%     generator         - instance generator
%     nSamples          - how large the sample to test on is
%     nItems            - how many items each instance has
%     policiesToRun     - a list of policies to test
%     exact             - compare to exact solutions found by binary programming
%     sampleFolder      - folder location of pregenerated instances
% Output:
%   valueFunctionEstimations        - a dictionary with results for each policy
%   valueFunctionTildeEstimations   - a dictionary with results in the approximated MDP used to generate policies
arguments
    policiesToRun
    generator 
    nSamples = 1000
    nItems = 200
    exact = 0
    sampleArray = []
end

%some preprocessing
if isempty(sampleArray)
    instances = cell(1,nSamples);
    for i = 1:nSamples
        instances{i} = instance('generator',generator,'nItems',nItems);
        while instances{i}.nItems == 0 %just a safety measure for weird instance generators.
            instances{i} = instance('generator',generator,'nItems',nItems);
        end
    end
else
    hold = load(sampleArray);
    instances = hold.savedKPinstances;
end
s=0; %make it a global variable
updateFrequency = 10;

if exact
    policiesToRunHold = [policiesToRun, "exact"];
    valueFunctionEstimations = dictionary(policiesToRunHold, struct);
else
    valueFunctionEstimations = dictionary(policiesToRun, struct);
end
valueFunctionTildeEstimations = dictionary;
if exact
    exactTic = tic;
    timeUpdate = updateFrequency;
    for instanceNumber = 1:nSamples
        getInstance
        sCopy = s;
        tic
        exactObj = sCopy.GurobiSolve;
        valueFunctionEstimations("exact").time(instanceNumber)      = toc;
        valueFunctionEstimations("exact").objective(instanceNumber) = exactObj;
        if toc(exactTic)>timeUpdate
            fprintf('Exact solution for %3.2f%% samples computed. This took %2.0f seconds \n', 100*instanceNumber/nSamples, toc(exactTic))
            timeUpdate = timeUpdate + updateFrequency;
        end
    end
    valueFunctionEstimations("exact").valueFunctionEstimate = sum(valueFunctionEstimations("exact").objective) / nSamples;
end
for policyCell = policiesToRun
    policy = policyCell{1};
    fprintf(['Running ', policy, ':\n'])
    nonGreedyCounter = 0;
    totalActions     = 0;
    sampledActions = 0;
    if strcmp(policy,'greedy')
        pi = @(s,t) max(1); %this is done such that asking for a second output does not raise an error
    elseif strcmp(policy, 'k-uniform')
        pi = @(s,t) max(randi(min(3,s.nItems),1));
    else
        pi = @(s,t) piPrime(s,t,policy);
        pi(s,1);%load the policy table already %assume an instance is loaded
        valueFunctionTildeEstimations{policy} = struct;
        load(policy, 'stateMappingTilde', 'valueFunctionPiTilde', 'valueFunctionPiTildePrime', 'abstractStates')
        valueFunctionPiTildePrime = valueFunctionPiTildePrime{1};
        objective_piTilde      = zeros(1,nSamples);
        objective_piTildePrime = zeros(1,nSamples);
    end
    times      = zeros(1,nSamples);
    objectives = zeros(1,nSamples);
    sampleds   = zeros(1,nSamples);
    policyTic = tic;
    timeUpdate = updateFrequency;
    for instanceNumber = 1:nSamples
        if toc(policyTic)>timeUpdate
            fprintf('%3.2f%% samples completed. This took %2.0f seconds \n', 100*instanceNumber/nSamples, toc(policyTic))
            timeUpdate = timeUpdate + updateFrequency;
        end
        getInstance
        sCopy = s;
        fullySampled = 1;
        if ~(strcmp(policy,'greedy') || strcmp(policy, 'k-uniform'))
            [~, stateIndex] = abstractStates.f(sCopy, 1);
            index = getQuickIndex(stateIndex, stateMappingTilde{1}, numel(stateMappingTilde{1}));
            if isempty(index) %Could happen if the new state is not encountered in the MonteCarlo Simulation
                objective_piTilde(instanceNumber) = NaN;
                objective_piTildePrime(instanceNumber) = NaN;
            else
                objective_piTilde(instanceNumber)      =  valueFunctionPiTilde(index);
                objective_piTildePrime(instanceNumber) =  valueFunctionPiTildePrime(index);
            end
        end
        tic
        for t = 1:nItems
            %piAction = pi(sCopy, t);
            [piAction, sampled] = pi(sCopy,t);
            if piAction ~=1
                nonGreedyCounter = nonGreedyCounter + 1;
            end
            totalActions = totalActions + 1;
            sampledActions = sampledActions + sampled;
            fullySampled = min(fullySampled, sampled);
            sCopy = sCopy.selectItem(piAction);
            if sCopy.nItems == 0
                times(instanceNumber)      = toc;
                objectives(instanceNumber) = sCopy.totalVal;
                sampleds(instanceNumber)   = fullySampled;
                break
            end
        end
    end
    valueFunctionEstimations(policy).time      = times;
    valueFunctionEstimations(policy).objective = objectives;
    valueFunctionEstimations(policy).sampled   = sampleds;
    valueFunctionEstimations(policy).valueFunctionEstimate = sum(valueFunctionEstimations(policy).objective) / nSamples;
    valueFunctionEstimations(policy).totalActions = totalActions;
    valueFunctionEstimations(policy).sampledActions = sampledActions;
    valueFunctionEstimations(policy).nonGreedyCounter = nonGreedyCounter;
    if ~(strcmp(policy,'greedy') || strcmp(policy, 'k-uniform'))
        valueFunctionTildeEstimations{policy}.piTilde      = objective_piTilde;
        valueFunctionTildeEstimations{policy}.piTildePrime = objective_piTildePrime;
    end
    if exact
        valueFunctionEstimations(policy).optGap = ( (valueFunctionEstimations(policy).valueFunctionEstimate / valueFunctionEstimations("exact").valueFunctionEstimate));
    end
    clear piPrime
end

    function getInstance
        s = instances{instanceNumber};
    end
    
    function [action,sampled] = piPrime_kuniform(s,t)
        [action, sampled] = piPrime(s,t,policy);
        if sampled==0
            action = randi(min(3,s.nItems),1);
        end
    end
end



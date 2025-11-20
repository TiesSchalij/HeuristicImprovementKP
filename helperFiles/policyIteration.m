function [piTildePrime, expectedPiTildePrimeValueFunction, valueFunctionPiTildePrime, varargout] = policyIteration(Pt, rsa, abstractStates, totalVisits, stateMapping, estimationCutoff, oldPolicyLocation)
% policyIteration performs policy iteration in an approximated TMDP. 
% Input:
%     Pt - transition dynamics
%     rsa - reward function
%     abstractStates - state space partition
%     totalVisits - how often each state was visited in the Monte Carlo simulation
%     stateMapping - list of state mapping
%     estimationCutoff - when is a state considered undersampled
%     oldPolicyLocation - where to find the old policy if this is a second round of improvement
% Output:
%     piTildePrime - an optimal policy in the approximated TMDP
%     expectedPiTildePrimeValueFunction - expected value function of pi tilde prime
%     valueFunctionPiTildePrime - value function of pi tilde prime, in accordance with statemapping
%     stateMapping - new state mapping (updated if this is the second round of improvement)
arguments
    Pt
    rsa
    abstractStates
    totalVisits
    stateMapping
    estimationCutoff = 0; %minimum amount of state visits to be able to deviate
    oldPolicyLocation = []; %only if second iteration of policy improving
end
T = abstractStates.T;
if ~(isempty(oldPolicyLocation) || strcmp(oldPolicyLocation, 'greedy'))
    if ~strcmp(oldPolicyLocation, 'k-uniform')
    oldPolicy = load(oldPolicyLocation);
    oldStateMapping = oldPolicy.stateMapping;
    oldPi = oldPolicy.piTildePrime;
    newStateMapping = cell(1,T);
    end
end


piTildePrime = cell(1,T);
valueFunctionPiTildePrime = cell(1,T);
valueFunctionPiTildePrime{T} = zeros(numel(rsa{T})/abstractStates.ks(T),1);
%last_t = find(cellfun('isempty', Pt),1, 'first') - 1
last_t = T;
valueFunctionPiTildePrime{last_t + 1} = 0; % terminal

for t = last_t:-1:1
    Qvalues_t = Pt{t} * valueFunctionPiTildePrime{t + 1} + rsa{t}; %calculate Qvalues
    nSA = numel(Qvalues_t);
    k = abstractStates.ks(t);
    Qvalues_t = reshape(Qvalues_t, k,nSA/k); %reshape from column vector to k by nSA/k matrix
    [valueFunction_t, policy] = max(Qvalues_t,[],1); %evaluate all Qvalues
    poorlyEstimatedStates = full(totalVisits{t})' < estimationCutoff;  %get states that are poorly sampled
    if ~(isempty(oldPolicyLocation) || strcmp(oldPolicyLocation, 'greedy'))
        if strcmp(oldPolicyLocation, 'k-uniform')
            wellEstimatedStates = full(totalVisits{t})' >= estimationCutoff;
            newStateMapping{t} = stateMapping{t}(wellEstimatedStates); %remove poorly estimated states, so they don't show up, and we can take a random action
            policy = policy(wellEstimatedStates);
        else
            wellEstimated = ~poorlyEstimatedStates;%don't change these states
            originalIndexWellEstimated = stateMapping{t}(logical([wellEstimated,0]));
            [statesToAddOriginalIndex,statesToAddNewIndex] = setdiff(oldStateMapping{t}(1:end-1),originalIndexWellEstimated); %original index of states that have to be added, and their new Index
            actionsForStatesToAdd=oldPi{t}(statesToAddNewIndex);%get action
            fprintf('In time step %d, we added %d states from the previous policy.\nOf which %d were non-greedy.\n\n', t, numel(statesToAddOriginalIndex), sum(actionsForStatesToAdd~=1))
            newStateMapping{t} = [stateMapping{t}(1:end-1), statesToAddOriginalIndex];
            policy = [policy, actionsForStatesToAdd'];
            newIndexPoorleyEstimatedStates = find(poorlyEstimatedStates);
            newIndexPoorleyEstimatedStates = setdiff(newIndexPoorleyEstimatedStates, statesToAddNewIndex);
            poorlyEstimatedStates = zeros(size(policy));
            poorlyEstimatedStates(newIndexPoorleyEstimatedStates) = 1;
        end
    end
    if ~strcmp(oldPolicyLocation, 'k-uniform')
        poorlyEstimatedStates = logical(poorlyEstimatedStates);
        valueFunctionPoorlyEstimated = Qvalues_t(1,poorlyEstimatedStates(1:nSA/k)); %don't deviate from greedy, so take first action
        policy(poorlyEstimatedStates) = 1;
        valueFunction_t(poorlyEstimatedStates(1:nSA/k)) = valueFunctionPoorlyEstimated;
    end
    valueFunction_t = [valueFunction_t,0];%#ok add terminal state

    valueFunctionPiTildePrime{t} = valueFunction_t';
    %policy and statemapping will mismatch from valuefunction and
    %transitions and rewards and everything
    if ~(isempty(oldPolicyLocation) || strcmp(oldPolicyLocation, 'greedy'))
        [newStateMappingHold,sortingIndex] = sort(newStateMapping{t});
        policy = policy(sortingIndex);
        newStateMapping{t} = newStateMappingHold;
    end
    piTildePrime{t} = policy';
end
expectedPiTildePrimeValueFunction = sum((valueFunctionPiTildePrime{1}(1:end-1,1)).*totalVisits{1} /sum(totalVisits{1})); %calculate the expected value function over all non terminal states
if ~(isempty(oldPolicyLocation) || strcmp(oldPolicyLocation, 'greedy'))
    varargout{1} = newStateMapping;
else
    varargout{1} = stateMapping;
end


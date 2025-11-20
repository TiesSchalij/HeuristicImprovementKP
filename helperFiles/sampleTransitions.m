function [Pt, rsa, stateMapping, totalVisits, expectedPiValueFunction, goesToTerminal, valueFunctionPiTilde] = sampleTransitions(abstractStates, policyName, nSamples, folderLocation, generator, nItems)
% sampleTransitions performs a Monte Carlo simulation to estimate the approximated TMDP
% Input:
%     abstractStates      - the state space partition
%     policyName          - the initial policy
%     nSamples            - how many runs of the Monte Carlo
%     folderLocation      - where to find a pregenerated sample for the Monte Carlo (optional)
%     generator           - what generator to use for the sample
%     nItems              - how many items in each instance
% Output:
%     Pt                        - the approximated transition dynamics
%     rsa                       - the approximated reward function
%     stateMapping              - a mapping from instance index to used instance index
%     totalVisits               - how often each approximated state has been sampled
%     expectedPiValueFunction   - the expected value function of pi 
%     goesToTerminal            - if a state can go to the terminal state
%     valueFunctionPiTilde      - the value function of pi tilde

arguments
    abstractStates
    policyName          =[]
    nSamples            =[]
    folderLocation      =[]
    generator           =[]
    nItems              =[]
end
%% Some preprocessing
tic
if isempty(policyName) || strcmp(policyName, 'greedy')%assume greedy policy
    policy = @(s,t,k) 1;
elseif strcmp(policyName, 'k-uniform')
    policy = @(s,t,k) randi(min(k,s.nItems),1);
else
    policy = @(s,t,k) piPrime(s,t,policyName);
end
if isempty(folderLocation)%none given --> use generator
    generateInstances = 1;
else
    generateInstances = 0;
    folder = dir(folderLocation);
    nFiles = numel(folder);
    currentArray = 3;%skip first 2 because of folder things
    arrayName = folder(currentArray).name;
    arrayWithInstances = open(strcat(folderLocation,'\', arrayName));
    arrayWithInstances = arrayWithInstances.savedKPinstances;
    nInstancesInArray = numel(arrayWithInstances);
    instanceNumberForFile = 1;
end
nSamples = min(5E6, nSamples); %5 million because I know we have around 3.6 or 4.5 million instances for the BPP testing
doneSampling  = 0;
instanceNumber = 0;
s=0;
expectedPiValueFunction = 0;
T = abstractStates.T;
last_t = 0;
timeUpdate = 10;%s
updateFrequency = 10;%s
%% Main loop
stateActionList = cell(1,T);
nextStateList   = cell(1,T);
actuallyUsedSA  = cell(1,T);
rewardSApairs   = cell(1,T);
goesToTerminal  = zeros(1,T);
sentFromAt_t    = cell(1,T);
for t = 1:T %preload storage
    stateActionList{t} = zeros(1,nSamples*abstractStates.ks(t)); %make a list to store all stateaction pairs at time t
    nextStateList{t}   = zeros(1,nSamples*abstractStates.ks(t)); %make a list to store all the next states
    actuallyUsedSA{t}  = zeros(1,nSamples*abstractStates.ks(t)); %store whether a SxA pair is valid
    rewardSApairs{t}   = zeros(1,nSamples*abstractStates.ks(t));
    sentFromAt_t{t}    = zeros(1,nSamples);
    initialStateInd    = zeros(1,nSamples);
    objectivesPi       = zeros(1,nSamples);
end


while ~doneSampling
    instanceNumber = instanceNumber + 1;
    if toc>timeUpdate
        fprintf('The Monte Carlo simulation is %.1f%% complete. This took %.1f seconds \n', instanceNumber*100/nSamples, toc)
        timeUpdate = timeUpdate + updateFrequency;
    end
    getInstance
    state = s;
    [oldState,oldStateInd] = abstractStates.f(state,1);
    initialStateInd(instanceNumber) = oldStateInd;
    for t = 1:T-1
        sentFromAt_t{t}(instanceNumber) = oldStateInd;
        k = abstractStates.ks(t);
        policyAction = policy(state,t,k);
        for action = 1:k
            if oldState(action) == abstractStates.nRegions(t) %no item here --> no valid action
                continue
            end
            sNew = state.selectItem(action); %get s'
            [newState, newStateIndex] = abstractStates.f(sNew, t+1); %and its index
            stateActionList{t}((instanceNumber-1)*k + action) = (oldStateInd-1)*k + action; %store unique value for state-action pair
            actuallyUsedSA{t}((instanceNumber-1)*k + action) = true;
            nextStateList{t}((instanceNumber-1)*k + action) = newStateIndex;
            rewardSApairs{t}((instanceNumber-1)*k + action) = state.values(action);
            if action == policyAction
                nextIterationStorageInstance = sNew;
                nextIterationStorageState = newState;
                nextIterationStorageIndex = newStateIndex;
            end
        end
        state = nextIterationStorageInstance;
        oldStateInd = nextIterationStorageIndex;
        oldState    = nextIterationStorageState;
        if state.nItems == 0%knapsack is full
            goesToTerminal(t) = 1;
            expectedPiValueFunction = expectedPiValueFunction + state.totalVal;
            objectivesPi(instanceNumber) = state.totalVal;
            last_t = max(last_t, t);
            break %no need to continue with the instance
        end % done with all action in one time period
    end % done with one instance

    if instanceNumber>=nSamples
        doneSampling = 1;
    end
end % done with all sampling
expectedPiValueFunction = expectedPiValueFunction / instanceNumber;
fprintf('The Monte Carlo simulation is 100 %% complete. This took %.1f seconds. In total %d samples were used \n\n', toc, instanceNumber)

%% Post processing
Pt = cell(1,T);
rsa= cell(1,T);
statesVisitedList = cell(1,T);
stateMapping = cell(1,T);
stateMapping2 = cell(1,T);
SAindexForP = cell(1,T);
nextStatesForP  = cell(1,T);
totalVisits = cell(1,T);

tic
for t = 1:last_t+1
    k = abstractStates.ks(t);
    statesVisitedList_temp = (stateActionList{t}(1:k:nSamples*k));
    statesVisitedList_temp(statesVisitedList_temp==0) = [];
    statesVisitedList{t} = ((statesVisitedList_temp - 1 )/k)+1; %convert to state index
    stateMapping{t} = unique(statesVisitedList{t});
    stateMapping{t}(stateMapping{t}==0) = [];

    stateMapping2{t} = unique(sentFromAt_t{t});
    stateMapping2{t}(stateMapping2{t}==0) = [];
    if ~isempty(stateMapping{t})
        if stateMapping{t}(end) < abstractStates.maxStates(t) %terminal state is not included
            %terminalStateAddedFor_t = t
            stateMapping{t} = [stateMapping{t}, abstractStates.maxStates(t)]; %include terminal state
        end
    else %still add terminal state
        stateMapping{t} = abstractStates.maxStates(t);
    end
    SAindexForP{t} = stateActionList{t}(1:instanceNumber*k);
    counterForMapping = 1;
    lengthOfStateMapping = numel(stateMapping{t});
    for i =SAindexForP{t}(1:k:end)
        if i==0 %SA not sampled
            counterForMapping = counterForMapping +k;
            continue
        end
        state = (i-1)/k+1;
        mappedIndex = getQuickIndex(state, stateMapping{t},lengthOfStateMapping);
        for a = 1:k
            SAindexForP{t}(counterForMapping) = (mappedIndex-1)*k + a;
            counterForMapping = counterForMapping + 1;
        end
    end
    fprintf('State Action Pairs for time period %d converted. This took %.1f seconds.\n', t, toc)
end
tic
for t = 1:last_t
    k = abstractStates.ks(t);
    nextStatesForP{t} = nextStateList{t}(1:instanceNumber*k);
    counterForMapping = 1;
    lengthOfStateMapping = numel(stateMapping{t+1});
    for i = nextStatesForP{t}
        if i==0 %SA not sampled
            counterForMapping = counterForMapping +1;
            continue
        end
        indexForNextState = getQuickIndex(i, stateMapping{t+1},lengthOfStateMapping);
        if isempty(indexForNextState) %the next state is not sampled --> maybe send to terminal state? that's a great idea!
            indexForNextState = lengthOfStateMapping; %terminal state
        end
        nextStatesForP{t}(counterForMapping) = indexForNextState;
        counterForMapping = counterForMapping + 1;
    end
    fprintf('Next states for time period %d converted. This took %.1f seconds.\n', t, toc)
end
tic
for t = 1:last_t
    k = abstractStates.ks(t);
    rewardsForP = rewardSApairs{t}(1:instanceNumber*k);
    Pt_temp = sparse(SAindexForP{t}(logical(actuallyUsedSA{t})), nextStatesForP{t}(logical(actuallyUsedSA{t})), 1);
    totalStateActionVisits = sum(Pt_temp,2);
    totalVisits{t} = totalStateActionVisits(1:k:end);
    fastTerm = 1./totalStateActionVisits;
    weightOfVisit = zeros(size(SAindexForP{t}(logical(actuallyUsedSA{t}))));
    counter = 1;
    for stateAction = SAindexForP{t}(logical(actuallyUsedSA{t}))
        weightOfVisit(counter) = fastTerm(stateAction);
        counter = counter + 1;
    end
    rsa{t} = sparse(SAindexForP{t}(logical(actuallyUsedSA{t})),1,rewardsForP(logical(actuallyUsedSA{t})).*weightOfVisit,(numel(stateMapping{t})-1)*k ,1);
    Pt{t}  = sparse(SAindexForP{t}(logical(actuallyUsedSA{t})), nextStatesForP{t}(logical(actuallyUsedSA{t})),weightOfVisit, (numel(stateMapping{t})-1)*k, numel(stateMapping{t+1}));
    fprintf('P and r for time period %d created. This took %.1f seconds.\n', t, toc)
end

valueFunctionPiTilde = zeros(numel(stateMapping{1}),1);
stateCounter = 1;
for state = stateMapping{1}
    samplesInState = initialStateInd==state;
    nSamplesInState = sum(samplesInState);
    valueFunctionPiTilde(stateCounter) = sum(objectivesPi(samplesInState))/nSamplesInState;
    stateCounter = stateCounter + 1;
end

%% Functions
    function getInstance
        if generateInstances
            s = instance('generator',generator, 'nItems',nItems);
        else %get next instance from folder
            s = arrayWithInstances{instanceNumberForFile};
            instanceNumberForFile = instanceNumberForFile +1;
            if instanceNumberForFile > nInstancesInArray %move to next array
                currentArray = currentArray + 1;
                if currentArray > nFiles % out of samples
                    doneSampling = 1;
                else
                    arrayName = folder(currentArray).name;
                    arrayWithInstances = open(strcat(folderLocation,'\', arrayName));
                    arrayWithInstances = arrayWithInstances.savedKPinstances;
                    nInstancesInArray = numel(arrayWithInstances);
                    instanceNumberForFile = 1;
                end
            end


        end

    end
if ~isempty(policyName)
    clear piPrime
end
end
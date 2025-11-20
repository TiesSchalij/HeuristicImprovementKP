clear; clc;%this is a file to create improved policies
%% User Input
T = 5;                                          %number of timeperiods
k = 7;                                          %number of actions to choose from
nItems = 200;                                   %number of items
nRows = 6; nCols = 6;                           %parameters for the Grid approximation
%if generating new instances                                         
generator = [];        %any generator or []
%if using presampled instances
instanceFolder = 'KPfromBPP_200items_200lb700ub1000C';                            %where to find the instances or []

policyName = 'greedy';                          %initial policy
nSamples = 1E6;                                 %how many instances to use in the Monte Carlo Simulation
estimationCutoff = 5;                           %undersampling parameter for the Monte Carlo Simulation
saveResults = 1;                                %if you want to save the results
if saveResults
    saveName = strcat('piTildePrime-',num2str(nRows),'_',num2str(nCols),'_',num2str(k),'-','BPP_20_70','-',num2str(nItems),'items.mat'); %under what name to save the results
end

%% Program 
grids = cellfun(@(x) {0:(1/nRows):((nRows - 1)/nRows),0:(1/nCols):((nCols-1)/nCols)},cell(1,T),'UniformOutput',false); %uniformly spaces the columns and rows
abstractStates = abstractStatesClass(T, grids, k*ones(1,T));
[Pt, rsa, stateMappingTilde, totalVisits, expectedPiValueFunction, goesToTerminal, valueFunctionPiTilde] = sampleTransitions(abstractStates,policyName, nSamples, instanceFolder, generator, nItems);

[piTildePrime, expectedPiTildePrimeValueFunction,valueFunctionPiTildePrime, stateMapping] = policyIteration(Pt, rsa, abstractStates, totalVisits, stateMappingTilde, estimationCutoff,policyName);
fprintf('Pi prime is %2.3f%% better than Pi.\n', (expectedPiTildePrimeValueFunction - expectedPiValueFunction)/expectedPiValueFunction * 100)
if saveResults
    save(saveName, 'piTildePrime', 'stateMapping', 'stateMappingTilde', 'abstractStates', 'totalVisits', 'estimationCutoff', 'valueFunctionPiTildePrime','expectedPiTildePrimeValueFunction', 'Pt','rsa' , 'valueFunctionPiTilde','-v7.3')
end


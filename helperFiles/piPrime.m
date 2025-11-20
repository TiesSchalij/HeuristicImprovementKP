function [action, varargout] = piPrime(instance, t, tableName)
% takes as input an instance, a time period, and a .mat file with relevant
% information. Creates persistent variable so that the .mat does not have
% to be loaded in each time.
arguments
    instance
    t
    tableName
end
persistent table
persistent mapping
persistent abstractS
persistent unsampledActions
if isempty(unsampledActions)
    unsampledActions = 0;
end
if isempty(table)
    temp = load(tableName);
    table = temp.piTildePrime;
    mapping = temp.stateMapping;
    abstractS = temp.abstractStates;
    % fprintf(strcat('max number of unique states for: ', tableName(19:25)))
    % max(cellfun(@numel, mapping))
end
[~, stateIndex] = abstractS.f(instance, t);
index = getQuickIndex(stateIndex, mapping{t}, numel(mapping{t}));
if isempty(index)%state is not sampled
    if strcmp(tableName(13:20), 'kuniform') %bit hacky but alright
        action = randi(min(3, instance.nItems));
    else 
        action = 1;%just take greedy action
    end
    if nargout >= 2
        varargout{1} = 0;
    end
    unsampledActions = unsampledActions + 1;
else
    action = table{t}(index);
    if nargout >= 2
        varargout{1} = 1;
    end
end
end
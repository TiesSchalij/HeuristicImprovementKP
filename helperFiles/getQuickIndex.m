function index = getQuickIndex(state, stateMapping, len)
%stateMapping is sorted in ascending order
%state is an element of stateMapping
upperbound = len;
if isempty(stateMapping)
    index = [];
    return
elseif stateMapping(upperbound) == state
    index = upperbound;
    return
end
lowerbound = 1;
index = floor((upperbound+lowerbound)/2);
found = stateMapping(index) == state;
secondaryCheck = 0;
while ~found
    if stateMapping(index)> state %decrease guess
        upperbound = index;
    else %increase guess
        lowerbound = index;
    end
    index = floor((upperbound+lowerbound)/2);
    if stateMapping(index) == state
        found = true;
    end
    if (upperbound-lowerbound) <= 1
        if secondaryCheck
            index = [];
            return
        else
            secondaryCheck = 1;
        end
    end
end
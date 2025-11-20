classdef abstractStatesClass
    properties
        T
        grids
        ks
        gridDim
        nRegions
        maxStates
        mergedStates
        dimNumbers
    end

    methods
        function obj = abstractStatesClass(T, grids, ks)
            % Input:
            % T: maximum time period
            % grids: 1xT cell array containing 1x2 cell arrays with the value and weight dividers
            % k: the number of items each time period considers
            % gridDim: 1xT cell array containing the dimensions of each grid
            % nRegions: number of regions in each grid
            % maxStates: maximum number of states for each grid
            obj.T = T;
            obj.grids = grids;
            obj.ks = ks;
            obj.gridDim = cellfun(@obj.getGridSize, grids);
            obj.nRegions = cellfun(@prod, obj.gridDim) + 1;%obj.nRows * obj.nCols; +1 for if an item is missing (null region)
            obj.maxStates = obj.nRegions.^obj.ks; %most of these states do not exist
            
            obj.dimNumbers = cell(1, T);
            for t = 1:T
                obj.dimNumbers{t} = obj.nRegions(t).^(0:obj.ks(t)-1);
            end

        end

        function varargout = f(obj, instance, t)
            %maps a Knapsack instance to an abstract state
            % Input:
            %    instance: the Knapsack Instance
            %    t: what time period
            % Output:
            %    state: the regions of the best k items
            %    stateInd (optional): the linear index of that state
            grid = obj.grids{t};
            nItemsToCheck = min(obj.ks(t), instance.nItems);%in case there are not enough items
            itemValues = instance.values(1:nItemsToCheck);
            itemWeights= instance.weights(1:nItemsToCheck);
            vInd = sum(itemValues' >grid{1},2);
            wInd = sum(itemWeights'>grid{2},2);

            nRows = obj.gridDim{t}(1);

            state = ((wInd-1)*nRows + vInd)';
            state = [state, obj.nRegions(t)*ones(1,obj.ks(t)-nItemsToCheck)];%add the null region
            varargout{1} = state;
            if nargout > 1
                varargout{2} = obj.state2ind(state, t);
            end
        end

        function  ind = state2ind(obj, state, t)
            %converts a state [r1, r2, ..., rk] to a linear index
            indices = state-1;
            %dimNumbers = (obj.nRegions(t)).^(0:(obj.ks(t)-1)); %manually calculate sub2ind
            ind = sum(obj.dimNumbers{t}.*indices)+1;
        end

        function state = ind2state(obj, ind, t)
            %converts a linear index of a state to a state [r1, r2, ..., rk]
            k = obj.ks(t);
            state = NaN(1,k);%preallocate
            nR = obj.nRegions(t);
            for i = k-1:-1:1
                state(i+1) = floor((ind-1)/(nR^(i)));%ind-1 to avoid the issue with integer rounding
                ind = ind - state(i+1)*nR^(i);
            end
            state = state +1;
            state(1) = ind; %remainder is the first region


        end

        function gridDim = getGridSize(~, grid)
            gridDim = {cellfun(@numel, grid)};
        end
    end
end
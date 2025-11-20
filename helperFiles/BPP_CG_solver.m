function [nColsGeneratedBP, BPTime, nColsGeneratedHeuristic, heuristicTime, nColsGeneratedTotal, totalTime, varargout] =BPP_CG_solver(n, C, w, solver1, solver2, saveKPs,verbose)
arguments
    n
    C
    w
    solver1
    solver2 = []
    saveKPs = 0
    verbose = 1
end

%initialize
totalTime = tic;
A = eye(n);
b = ones(n,1);
c = ones(1,n);
prevSol = ones(n,1);
nColsGeneratedHeuristic = 0;
nColsGeneratedBP        = 0;
nColsGeneratedTotal     = -1;
heuristicTime  = 0;
BPTime         = 0;
newColumnFound = 1;
if ~isempty(solver2)
    totalXX = 0;totalXY = 0;totalYX = 0;totalYY = 0;
end

if saveKPs
    saveFolderKP        = 'KPfromBPP_200items_200lb700ub1000C';
    savedKPinstances = cell(1);
end


while newColumnFound
    nColsGeneratedTotal = nColsGeneratedTotal +1;
    newColumnFound = 0;
    %% MasterProblem

    BPPmodel.A= sparse(A);
    BPPmodel.obj = c;
    BPPmodel.rhs = b;
    BPPmodel.sense ='=';
    BPPmodel.modelsense = 'min';
    BPPparams.OutputFlag = 0;
    BPPmodel.start = prevSol;
    output = gurobi(BPPmodel, BPPparams);

    %% Subproblem
    v = output.pi;
    s = instance('v' , v', 'w', w/C);
    scaleFactor = max(v);
    if saveKPs
        savedKPinstances{nColsGeneratedTotal + 1} = s;
    end
    heuristicStart = tic;
    [subSol, XX, XY, YX, YY] = subroutine(s, solver1, solver2,scaleFactor);
    heuristicTime = heuristicTime + toc(heuristicStart);
    if ~isempty(solver2)
        totalXX = totalXX + XX;totalXY = totalXY + XY;totalYX = totalYX + YX;totalYY = totalYY + YY;
    end
    newColumn = zeros(n,1);
    newColumn(subSol) = 1;
    if round(1 - sum(v.*newColumn),12) < 0 %newColumn has negative reduced cost
        A = [A, newColumn];
        nColsGeneratedHeuristic = nColsGeneratedHeuristic + 1;
        c = [c,1];
        newColumnFound = 1;
    else %check with Gurobi
        BPstart = tic;
        [~,BPSol] = s.GurobiSolve(0,inf,0, subSol);
        BPTime = BPTime + toc(BPstart);
        newColumn = zeros(n,1);
        newColumn(BPSol) = 1;
        if round(1 - sum(v.*newColumn),12) < 0 %newColumn has negative reduced cost
            A = [A, newColumn];
            nColsGeneratedBP = nColsGeneratedBP + 1;
            c = [c,1];
            newColumnFound = 1;
        end
    end
    prevSol = [output.x; 0]; %add 0 for the new column
end
totalTime = toc(totalTime);
clear piPrime %just to be sure
if saveKPs
    arrayNumber = numel(dir(saveFolderKP)) -1;
    save(strcat(saveFolderKP,'\','instanceArray',string(arrayNumber) ),'savedKPinstances')
    fprintf('Saving instances successful\nNow %d arrays with KP instances saved\n', arrayNumber)
end
if verbose
    fprintf('The total runtime was: %.2f seconds, the exact solver took: %.2f seconds, and the subproblem solver took: %.2f seconds\n', totalTime, BPTime, heuristicTime)
    fprintf('Optimal solution found: %.6f bins needed\nThere were %d columns generated.\nThe subproblem was solved to optimality %d times.\n', output.objval, nColsGeneratedTotal, nColsGeneratedBP)
end

if ~isempty(solver2)
    varargout{1} = totalXX;varargout{2} = totalXY;varargout{3} = totalYX;varargout{4} = totalYY;
end

    function [sol, XX, XY, YX, YY] = subroutine(s,solver1,solver2,scaleFactor)
        r = round(1/scaleFactor,12);
        [sol1, obj1] = solver1(s);
        if isempty(solver2)
            sol = sol1;
            obj2 = 0;
        else
            [sol2, obj2] = solver2(s);
            if obj2>obj1
                sol = sol2;
            else
                sol = sol1;
            end
        end
        obj1 = round(obj1,12);obj2=round(obj2,12);
        XX = (obj1 <= r) & (obj2 <= r); XY = (obj1 <= r) & (obj2 > r);YX = (obj1 > r) & (obj2 <= r); YY = (obj1 > r) & (obj2 > r);
    end
end

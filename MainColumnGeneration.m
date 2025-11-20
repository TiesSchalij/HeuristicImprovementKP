clear; clc;
%% User Input:
nRuns = 100;                    %number of Bin Packing Problems to solve
nItems = 200;                   %number of items in each problem
BPParray = [];                  %filename for presampled BPP instances
saveKPs = 0;
LaTeXTable = 1;
saveBPPs = 0;

subroutines =["greedy",...
              "Improved",...
              "greedy+Improved"];

%% Program
if isempty(BPParray)
    BPPinstances = cell(1);
    for i = 1:nRuns
        [n, C, w] = BPPrandomInstance(1000, nItems);
        BPPinstance.n = n;
        BPPinstance.C = C;
        BPPinstance.w = w;
        BPPinstances{i} = BPPinstance;
    end
else
    load(BPParray)
end
results = dictionary(subroutines, struct);
for run = 1:nRuns
    fprintf('\nCommencing run number %d:\n', run)
    n = BPPinstances{run}.n;
    C = BPPinstances{run}.C;
    w = BPPinstances{run}.w;
    for subroutine = subroutines
        fprintf(append('\nRunning ',subroutine,' policy:\n'))
        if strcmp(subroutine, 'greedy')
            [nColsGeneratedBP, BPTime, nColsGeneratedHeuristic, heuristicTime, nColsGeneratedTotal, totalTime] = BPP_CG_solver(n, C, w, @greedy,[],0,1);
        elseif strcmp(subroutine, 'Improved')
            [nColsGeneratedBP, BPTime, nColsGeneratedHeuristic, heuristicTime, nColsGeneratedTotal, totalTime] = BPP_CG_solver(n, C, w, @improved,[],0,1);
        elseif strcmp(subroutine, 'greedy+Improved')
            [nColsGeneratedBP, BPTime, nColsGeneratedHeuristic, heuristicTime, nColsGeneratedTotal, totalTime,XX,XY,YX,YY] = BPP_CG_solver(n, C, w, @greedy,@improved,0,1);
        elseif strcmp(subroutine, 'greedy+ImprovedFailed')
            [nColsGeneratedBP, BPTime, nColsGeneratedHeuristic, heuristicTime, nColsGeneratedTotal, totalTime,XX,XY,YX,YY] = BPP_CG_solver(n, C, w, @greedy,@improvedFailed,0,1);
        end
        results(subroutine).nColsGeneratedBP(run) = nColsGeneratedBP;
        results(subroutine).BPTime(run) = BPTime;
        results(subroutine).nColsGeneratedHeuristic(run) = nColsGeneratedHeuristic;
        results(subroutine).heuristicTime(run) = heuristicTime;
        results(subroutine).nColsGeneratedTotal(run) = nColsGeneratedTotal;
        results(subroutine).totalTime(run) = totalTime;
        if strcmp(subroutine, 'greedy+Improved') | strcmp(subroutine, 'greedy+ImprovedFailed')
            results(subroutine).XX(run) = XX;results(subroutine).XY(run) = XY;results(subroutine).YX(run) = YX;results(subroutine).YY(run) = YY;
        end

    end
end
if saveBPPs
    save('BPPinstancesArray.mat', 'BPPinstances')
end
%%
for key = subroutines
    R = results(key);

    mBP = mean(R.nColsGeneratedBP);
    mH  = mean(R.nColsGeneratedHeuristic);
    mT  = mean(R.nColsGeneratedTotal);

    tBP = mean(R.BPTime);
    tH  = mean(R.heuristicTime);
    tT  = mean(R.totalTime);

    fprintf("\nSummary for %s:\n", (key));
    fprintf("Average total runtime: %.2f seconds\n", tT);
    fprintf("  Gurobi time: %.2f seconds\n", tBP);
    fprintf("  Subroutine time: %.2f seconds\n", tH);

    fprintf("Average number of Columns generated: %.2f total\n", mT);
    fprintf("  via Gurobi: %.2f\n", mBP);
    fprintf("  via Subroutine: %.2f\n\n", mH);
    if strcmp(key, 'greedy+Improved') | strcmp(key, 'greedy+ImprovedFailed')
        fprintf("Did the subroutines find a negative reduced cost column?\n")
        fprintf("   Greedy no, improved no: %.2f total\n", mean(R.XX))
        fprintf("   Greedy no, improved yes: %.2f total\n", mean(R.XY))
        fprintf("   Greedy yes, improved no: %.2f total\n", mean(R.YX))
        fprintf("   Greedy yes, improved yes: %.2f total\n\n", mean(R.YY))
    end
end

if LaTeXTable
    fprintf("\\begin{tabular}{lcccccc}\n");
    fprintf("\\toprule\n");
    fprintf("Subroutine & \\\\multicolumn{3}{c}{Columns Generated} & \\\\multicolumn{3}{c}{Time (s)} \\\\\n");
    fprintf("\\\\cmidrule(lr){2-4} \\\\cmidrule(lr){5-7}\n");
    fprintf(" & via Gurobi & via Subroutine & Total & Gurobi & Subroutine & Total \\\\\n");
    fprintf("\\midrule\n");

    for key = subroutines
        R = results(key);

        mBP = mean(R.nColsGeneratedBP);
        mH  = mean(R.nColsGeneratedHeuristic);
        mT  = mean(R.nColsGeneratedTotal);

        tBP = mean(R.BPTime);
        tH  = mean(R.heuristicTime);
        tT  = mean(R.totalTime);

        fprintf("%s & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f \\\\\n", ...
            key, mBP, mH, mT, tBP, tH, tT);
    end

    fprintf("\\bottomrule\n");
    fprintf("\\end{tabular}\n");
end

function [sol, obj] = greedy(s)
t=1;
sol = [];
obj = 0;
while s.nItems>0
    action = 1;
    sol = [sol, s.originalIndex(action)];
    obj = obj + s.values(action);
    s = s.selectItem(action);
    t=t+1;
end
end

function [sol, obj] = improved(s)
t=1;
sol = [];
obj = 0;
while s.nItems>0
    action = piPrime(s,t,'piTildePrime-12_12_3-BPP_20_70-200items.mat');
    if action > s.nItems
        break % negative value item
    else
        sol = [sol, s.originalIndex(action)];
        obj = obj + s.values(action);
        s = s.selectItem(action);
        t=t+1;
    end
end
end

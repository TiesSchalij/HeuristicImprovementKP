# HeuristicImprovementKP_Code
This repository contains the code used for the paper titled 'From Generic to Structure-Aware Heuristics through Markov Decision Processes'.  
All files are Matlab files. There is a dependency on Gurobi, this can be avoided, which will be explained later.  

There are 3 main files:  
MainPolicyCreation.m      |  Approximates a TMDP and performs policy iteration to obtain an improved heuristic.  
MainPerformanceTesting.m  |  Evaluates the performance of heuristics.  
MainColumnGeneration.m    |  Solves BPP relaxations using column generation using different heuristics for the subproblem.  

The rest of the files are helper files. Please ensure Matlab has access to these files when running the main files.  
All three main files should be self-explanatory, but here follows a brief instruction.  

## HOW TO CREATE A POLICY IN `MainPolicyCreation.m`:
- Specify T, N, and k that define the TMDP encompassing the KnapsackProblem (section 5.2).  
- Specify the number of rows and columns for the (m,n)-Grid-Approximation (section 5.3).  
- Specify what generator to use as D₀ (can be any of the Pisinger generators in the helperFiles) or give a folder location with presampled instances.  
- Choose an initial policy. This can be 'greedy', 'k-uniform', or any previously trained policy.  
- Specify parameters of the Monte Carlo simulation (section 5.4).  
- Finally, specify if and how the obtained policy should be saved.  

The Monte Carlo simulation prints progress updates. After Policy Iteration is complete the relative improvement is printed.  

## HOW TO EVALUATE A POLICY WITH `MainPerformanceTesting.m`:
- Specify what generator to use or where to find pregenerated instances.  
- Specify how many instances are tested and how many items they have.  
- Specify the heuristics you wish to test on the same instances.  
- Specify whether you wish to compare to the exact solution found via Gurobi (may be slow).  
  In case Gurobi is unavailable, in line 52 of `performanceEstimation.m` change  
  `sCopy.GurobiSolve` to `sCopy.BPsolve`. This uses the Matlab built-in `intlinprog` instead.  

The file prints summarizing statistics after it is finished.  
The final three lines are only valid if five policies are tested at once, like in section 5.4.2. Simply comment these out if this is not relevant.  
The trained policies for the Pisinger generators are given in the folder `policies`.  
The data used to obtain the results in the paper are given in the folder `testingData`.  
Using a policy on its corresponding dataset will yield the results presented in the paper.  

## HOW TO SOLVE A BPP RELAXATION WITH `MainColumnGeneration.m`:
- Specify how many BPP instances you wish to solve and how many items each instance has.  
- Specify where to find presampled BPP instances; if left empty, new instances will be generated.  
- Specify whether you want to save all intermediately generated KP instances.  
  If so, a folder under the name `KPfromBPP_200items_200lb700ub1000C` should be accessible  
  (this name can be changed in line 29 of `BPP_CG_solver.m`).  
- Specify whether you want to print a LaTeX table.  
- Specify whether you want to save the used BPP instances.  
- Specify which heuristics you wish to use for the subroutine.  
  The heuristic for "improved" must match the filename as in line 131.  

Summarizing statistics are printed at the end of the program.  

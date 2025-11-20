function [v,w] = Pf_subsetSum(nItems, R)
arguments
    nItems = 200
    R      = 10000
end
w = randi(R, 1, nItems);
v = w;
w = w/R;
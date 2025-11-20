function [v,w] = Pd_inverseStronglyCorrelated(nItems, R)
arguments
    nItems = 200
    R      = 10000
end
v = randi(R, 1, nItems);
w = v + R/10;
w = w/R;
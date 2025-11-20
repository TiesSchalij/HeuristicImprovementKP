function [v,w] = Pc_stronglyCorrelated(nItems, R)
arguments
    nItems = 200
    R      = 10000
end
w = randi(R,1,nItems);
v = w + R/10;
w = w/R;
function [v,w] = Pe_almostStronglyCorrelated(nItems, R)
arguments
    nItems = 200
    R      = 10000
end
w = randi(R, 1, nItems);
v_noise = randi([-R,R]/500,1,nItems) + R/10;
v = w + v_noise;
w = w/R;
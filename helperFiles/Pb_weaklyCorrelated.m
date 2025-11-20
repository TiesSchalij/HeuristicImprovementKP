function [v,w] = Pb_weaklyCorrelated(nItems, R)
arguments
    nItems = 200
    R      = 10000
end
w = randi(R, 1, nItems);
v_noise = randi([-R, R]/10, 1, nItems);
v = max(w + v_noise, 1);
w = w/R;
function [v,w] = Pa_uncorrelated(nItems, R)
arguments
    nItems = 200
    R      = 10000
end
hold = randi(R, 2, nItems);
v = hold(1,:);
w = hold(2,:);
w = w/R;
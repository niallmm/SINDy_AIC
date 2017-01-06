function [value, isterminal, direction] = events_lorenz(t,y);

value = abs(max(y))-1e5;
isterminal = 1;
direction = 0;


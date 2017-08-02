function [value, isterminal, direction] = events_diverge(t,y);

value = max(abs(y))-1e9;
isterminal = 1;
direction = [];


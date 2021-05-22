function y = smooth_sign(x, scale)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
% y = x;
% y(x>=0) = 1 - exp(-scale*x(x>=0));
% y(x<0) = exp(scale*x(x<0)) - 1;

y = tanh(x*scale);
end


function output = Saturate(input, lower, upper)
%SATURATE Saturates the input subject to bounds.
%   SATURATE(INPUT, LOWER, UPPER) is a vectorized function that, for each
%   element of the vector INPUT returns LOWER if INPUT < LOWER, UPPER if
%   INPUT > UPPER, and INPUT otherwise.

output = max(input, lower);
output = min(output, upper);

end


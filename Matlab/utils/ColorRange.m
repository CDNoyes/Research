% COLORRANGE Returns a vector C with NPOINTS points to be use in plotting.
% The range is optimized to change smoothly from COLOR1 to COLOR2

function c = ColorRange(colorA, colorB, nPoints)
for i=1:3
    c(1:nPoints,i) = linspace(colorA(i), colorB(i), nPoints);
end
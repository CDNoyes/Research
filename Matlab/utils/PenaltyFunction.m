%PENALTYFUNCTION Defines a penalty function to be used for soft constraints
%in optimal constrol problems.
%   PENALTYFUNCTION(U,U_LOW,U_HIGH,B) Implements an inverse tangent smooth
%   function that penalizes values of U outside the bounds [U_LOW,U_HIGH].
%   B is a shaping parameter that determines the output value for values of
%   U outside the bounds.

%   See "Entry Trajectory Planner for Higher Elevation Mars Landing" by
%   Benito et al. for reference.
%
%   Could potentially return a function instead of simply evaluating at U.
%   Although this is easily accomplished via anonymous functions, i.e.
%   PF = @(x) PenaltyFunction(x, -10,10,100)

function L = PenaltyFunction(u, u_low, u_high, b)

L = (atan(b*(-u+u_high))+atan(b*(u-u_low)))/sign(u_high);

end


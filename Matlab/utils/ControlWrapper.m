%CONTROLWRAPPER Helpful utility that allows control inputs to be either a
%function or a value (or vector of values).
%   CONTROLWRAPPER(TIME,CONTROL) outputs the result of the function
%   CONTROL evaluated at TIME if CONTROL is a function and simply outputs CONTROL
%   otherwise. This allows one to define polymorphic system dynamics - in
%   one case the control is a function of time (as in the case of
%   interpolating a discrete control parametrization) for use in
%   integration; in the second case, the control is simply a value of
%   vector of values for use in numerically determining derivatives.

function u = ControlWrapper(time,control)

if isa(control,'function_handle')
    u = control(time);
else
    u = control;

end
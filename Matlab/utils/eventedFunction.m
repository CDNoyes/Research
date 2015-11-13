% EVENTEDFUNCTION Takes system dynamics and trivializes them (sets all
% derivatives to 0) once a given condition is met.
%   MATLAB's event system can be laborious to use.
%   EVENTEDFUNCTION(T,X,FUN,CONDITION) simplifies the process by returning
%   the original system's derivatives if the condition has not been met,
%   and all zeros if the condition has been met.
%
%   This might not be worthwhile. It may be better to simply put the
%   condition in the actual derivative file. It does separate out the
%   conditional event from the actual simulation code so maybe it does
%   serve a purpose after all.

function dx = eventedFunction(t,x,fun,condition) 

if eval(condition); 
    dx = zeros(size(x)); 
else
    dx = fun(t,x);
end

end
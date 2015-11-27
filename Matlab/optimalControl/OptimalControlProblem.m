%OPTIMALCONTROLPROBLEM Creates a consistent structure to be utilized by
%methods that solve optimal control problems.
%
%The OCP structure includes the dynamics, lagrange and mayer cost
%functions, constraints on state and control, as well as bounds on time,
%state, and control that can be additionally be used to specify fixed/free
%initial/final states. Any fields that are unused should be input as empty
%brackets.

function ocp = OptimalControlProblem(ode, lagrange, mayer, constraints, bounds,hessians)

if nargin < 6 || isempty(hessians)
    ocp.hessians = [];
elseif ~isstruct(hessians)
    error('When provided, HESSIANS must be a structure with at least one of the following fields: [lagrange,mayer,dynamics]')
else
    ocp.hessians = hessians; %Should check for field names, add them if not there.
end

ocp.dynamics = ode;
ocp.cost.lagrange = lagrange;
ocp.cost.mayer = mayer;
ocp.constraints = constraints; %Fine if empty, each individual method decides how to proceed.

%Example of the bounds structure:
% bounds.lower.control = -1;
% bounds.upper.control = 1;
% bounds.upper.initialState = [x0;30];
% bounds.lower.initialState = [x0;-30];
% bounds.upper.state = [10,10,30];
% bounds.lower.state = [-10,-10,-30];
% bounds.upper.finalState = [xf;30];
% bounds.lower.finalState = [xf;-30];
% bounds.upper.finalTime = tf;
% bounds.lower.finalTime = tf;

% We can get a lot of information from the bounds structure
ocp.bounds = bounds;
if (bounds.upper.finalTime-bounds.lower.finalTime) == 0
    ocp.fixedFinalTime = true;
else
    ocp.fixedFinalTime = false;
end

ocp.dimension.state = length(bounds.upper.initialState);
ocp.dimension.control = length(bounds.upper.control);
ocp.dimension.total = ocp.dimension.state+ocp.dimension.control;
%If an initial or final state of an element in the state has equal upper
%and lower bounds, then it should included as an equality constraint.
%Otherwise, it is included as two inequality constraints.
ocp = findFreeParameters(ocp);

% Not all formulations will use the adjoint but for those that do:
ocp.dimension.adjoint = ocp.dimension.state-length(~ocp.free.finalState);
end

function OCP = findFreeParameters(OCP)
%Based on bounds, determine which initial and final conditions are free and
%should therefore be included in the optimization vector.
% OCP.free.initialState = find(OCP.bounds.upper.initialState-OCP.bounds.lower.initialState);
OCP.free.initialState = (OCP.bounds.upper.initialState~=OCP.bounds.lower.initialState);
OCP.free.initialStateBool = ~isempty(OCP.free.initialState);
% OCP.free.finalState = find(OCP.bounds.upper.finalState-OCP.bounds.lower.finalState);
OCP.free.finalState = (OCP.bounds.upper.finalState~=OCP.bounds.lower.finalState);
OCP.free.finalStateBool = ~isempty(OCP.free.finalState);

end
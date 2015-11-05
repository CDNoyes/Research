%OPTIMALCONTROLPROBLEM Creates a consistent structure to be utilized by
%methods that solve optimal control problems.
%
%The OCP structure includes the dynamics, lagrange and mayer cost
%functions, constraints on state and control, as well as bounds on time,
%state, and control that can be additionally be used to specify fixed/free
%initial/final states. Any fields that are unused should be input as empty
%brackets.

function ocp = OptimalControlProblem(ode, lagrange, mayer, constraints, bounds)

ocp.dynamics = ode;

% if isempty(lagrange)
%     ocp.cost.lagrange = @(x) 0;
% else
ocp.cost.lagrange = lagrange;
% end

% if isempty(mayer)
%     ocp.cost.mayer = @(x) 0;
% else
ocp.cost.mayer = mayer;
% end

ocp.constraints = constraints; %Fine if empty, each individual method decides how to proceed.

% problem.bounds.lower.control = -1;
% problem.bounds.upper.control = 1;
% problem.bounds.upper.initialState = [x0;30];
% problem.bounds.lower.initialState = [x0;-30];
% problem.bounds.upper.state = [10,10,30];
% problem.bounds.lower.state = [-10,-10,-30];
% problem.bounds.upper.finalState = [xf;30];
% problem.bounds.lower.finalState = [xf;-30];
% problem.bounds.upper.finalTime = tf;
% problem.bounds.lower.finalTime = tf;

% We can get a lot of information from the bounds structure
ocp.bounds = bounds;
if (bounds.upper.finalTime-bounds.lower.finalTime) == 0
    ocp.fixedFinalTime = true;
else
    ocp.fixedFinalTime = false;
end

ocp.dimension.state = length(bounds.upper.initialState);
ocp.dimension.control = length(bounds.upper.control);

%If an initial or final state of an element in the state has equal upper
%and lower bounds, then it should included as an equality constraint.
%Otherwise, it is included as two inequality constraints.
ocp = findFreeParameters(ocp);

% Not all formulations will use the adjoint but for those that do:
ocp.dimension.adjoint = length(ocp.fixed.finalState);
end

function OCP = findFreeParameters(OCP)
%Based on bounds, determine which initial and final conditions are free and
%should therefore be included in the optimization vector.
OCP.free.initialState = find(OCP.bounds.upper.initialState-OCP.bounds.lower.initialState);
OCP.free.initialStateBool = ~isempty(OCP.free.initialState);
OCP.free.finalState = find(OCP.bounds.upper.finalState-OCP.bounds.lower.finalState);
OCP.free.finalStateBool = ~isempty(OCP.free.finalState);

end
%PLANNERDYNAMICS Integrable dynamics for use in the high elevation planner.
%   PLANNERDYNAMICS(TIME,STATE,BANKANGLE,PLANETMODEL,VEHICLEMODEL) is a
%   Matlab ode-solver compliant function that returns trivial dynamics when
%   the parachute deployment constraints on altitude and velocity have been
%   met and standard entry dynamics otherwise.

function dX = PlannerDynamics(t,x,sigma,planetModel,vehicleModel)

r_eq = planetModel.radiusEquatorial;

%Check parachute constraints - stop if slow enough or too low.
hmin = 6; %km
vmax = 420; %m/s
if x(4) < vmax || (x(1)-r_eq)/1000 < hmin
    dX = zeros(size(x));
else
    [g,L,D] = EntryForces(x,planetModel,vehicleModel);
    dX = EntryDynamics(x,sigma,g,L,D);
end

end
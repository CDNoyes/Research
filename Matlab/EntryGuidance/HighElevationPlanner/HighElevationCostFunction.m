% HIGHELEVATIONCOSTFUNCTION Objective function used in optimization to
% compute entry trajectories with high elevations.
%   HIGHELEVATIONCOSTFUNCTION(P,PLANETMODEL,VEHICLEMODEL) computes the cost
%   of a trajectory corresponding to a bank angle profile parametrized by
%   the three switching times in P.

function [J,t,x] = HighElevationCostFunction(p,planetModel,vehicleModel,DR,CR)

%Drive the optimization away from negative switching times and unordered
%times
J = (checkFeasibility(p));
if J > 1e3
    [t,x(1:6)] = deal(0);
    return
end

%Optimal control results show that the switching times are never outside
% [50,220] so we use this info to reduce the search space
if any(p<50) || any(p>220)
    J = 2e5*sum(abs(p-Saturate(p,50,220)));
    [t,x(1:6)] = deal(0);
    return
end
  
t1 = p(1);
t2 = p(2);
t3 = p(3);

dtr = pi/180;

x0 = [3540e3; -90.07*dtr; -43.90*dtr; 5505; -14.15*dtr; 4.99*dtr];

tf = 350; %Just needs to be long enough

r_eq = planetModel.radiusEquatorial;      % equatorial radius, m

sigma_min = 18.19*dtr;
sigma_max = 87.13*dtr;
fun = @(t) BankAngleProfile(t,t1,t2,t3,sigma_min,sigma_max);

opt = odeset('RelTol',1e-7,'AbsTol',1e-7);
[t,x] = ode45(@(T,X) PlannerDynamics(T,X,fun(T),planetModel,vehicleModel),linspace(0, tf,1000), x0,opt);

%Cost function weights:
k_h = 1e-7;
k_gamma = 0;
k_d = 1;


h = x(end,1) - r_eq;
phi = x(end,3);
theta = x(end,2);
gamma = x(end,5);


%Distance metric using range:
[dr,cr] = Range(x0(2),x0(3),x0(6),theta,phi);
d = norm([dr-DR,cr-CR]);
J = (-h*k_h + k_gamma*gamma^2 + k_d*d);

end

function cost = checkFeasibility(t)
%Avoid computing the actual cost if the chosen times are infeasible. The
%cost will be 0 if the three switching times are feasible.
    sig = (-diff(t));       
    cost = 1e6*(sum(sig+abs(sig))+sum(-t+abs(t)));
    if cost %#ok<BDLGI>
        cost = max(cost,1e7);
    end

end
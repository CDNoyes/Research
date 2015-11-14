function dX = PlannerDynamics(t,x,sigma,x0,DR)

% [currentDR,CR] = Range(x0(2),x0(3),x0(6),x(2),x(3));

%Compute range and compare to DR
%Also check parachute constraints - stop if too slow or too low.
r_eq = 3397e3;
hmin = 6; %km
vmin = 480; %m/s
if x(4) < vmin || (x(1)-r_eq)/1000 < hmin %|| currentDR > DR
    dX = zeros(size(x));
else
    [g,L,D] = EntryForces(x);
    dX = EntryDynamics(x,sigma,g,L,D);
end

end
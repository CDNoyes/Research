% TRAJLENGTH Computes the range (in radians) on the surface of a spherical
% body along the great circle connecting the initial and final locations.

function s = TrajLength(trajsum)
time = trajsum.time;
x = trajsum.state;
sdot = @(t,s) interp1(time,x(:,4).*cos(x(:,5))./x(:,1),t);

[t,s] = ode45(sdot,time,0);

end
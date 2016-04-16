function s = TrajLength(trajsum)
time = trajsum.time;
x = trajsum.state;
sdot = @(t,s) interp1(time,x(:,4).*cos(x(:,5)),t);

[t,s] = ode45(sdot,time,0);



end
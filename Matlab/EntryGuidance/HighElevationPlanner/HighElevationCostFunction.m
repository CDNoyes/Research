function [J,t,x] = HighElevationCostFunction(p)


t1 = p(1);
t2 = p(2);
t3 = p(3);

tf = 300;

r_eq = 3397e3;      % equatorial radius, m

dtr = pi/180;
sigma_min = 18.19*dtr;
sigma_max = 87.13*dtr;
fun = @(t) BankAngleProfile(t,t1,t2,t3,sigma_min,sigma_max);

dyn = @(t,x) EntryDynamics(x, fun(t));

x0 = [-90.07*dtr; -43.90*dtr; 3540e3; 4.99*dtr; -14.15*dtr; 5505];
opt = odeset(opt,'RelTol',1e-8,'AbsTol',1e-8);
[t,x] = ode45(dyn,[0 tf], x0,opt,fun);

k_h = 5;
k_gamma = (0.1*dtr)^2;
k_d = 1;


h = x(:,3) - r_eq;
phi = x(:,2);
theta = x(:,1);
DR = 750;
CR = 0;
[theta_T,phi_T] = FinalLatLon(x0(1),x0(2),x0(4),DR,CR);
d = 2*r_eq*asin(sqrt( sin(0.5*(phi_T-phi)).^2 + cos(phi_T).*cos(phi).*sin(0.5*(theta_T-theta)).^2 ));

J = -h(end)*kh + k_gamma*gamma^2 + k_d*d;

end
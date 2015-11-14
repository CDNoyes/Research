%ENTRYDYNAMICS Computes the state derivatives of an entry vehicle.
%   ENTRYDYNAMICS(X,SIGMA, SCALEFACTOR)

function dX = EntryDynamics(x,sigma,g,L,D)

% omega_p = 7.095e-5; % angular rate of planet rotation, rad/s

%State Variables:
% x = [ r theta phi V gamma psi ]

% theta = x(2);   %longitude, rad
phi = x(3);     %latitude, rad
r = x(1);       %vehicle radius, m
psi = x(6);     %heading angle, rad, 0 -> due East
gamma = x(5);   %flight path angle, rad
V = x(4);       %velocity, m/s

r_dot = V*sin(gamma);
V_dot = -D-g*sin(gamma);
theta_dot = V/r*cos(gamma)*cos(psi)/cos(phi);
gamma_dot = (L*cos(sigma) - (g-V^2/r)*cos(gamma))/V;% + 2*omega_p*cos(psi)*cos(phi);
phi_dot = V/r*cos(gamma)*sin(psi);
psi_dot = 1/(V*cos(gamma))*(L*sin(sigma)-V^2/r*cos(gamma)^2*cos(psi)*tan(phi));% + 2*omega_p*(tan(gamma)*sin(psi)*cos(phi)-sin(phi));

dX = [r_dot; theta_dot; phi_dot;  V_dot; gamma_dot; psi_dot ];
end
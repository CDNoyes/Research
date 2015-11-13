%Generic Entry dynamics for use in reference tracking problems where the
%equations need to be called twice, once for the reference and once for the
%perturbed/flown trajectory.

function [dX,D,hs] = EntryDynamics(x,sigma,ScaleFactor,Perturbation)
%Constants:
r_eq = 3397e3;      % equatorial radius, m
% omega_p = 7.095e-5; % angular rate of planet rotation, rad/s
mu = 4.2830e13;     % gravitational parameter, m^3/s^2
if nargin < 3 || isempty(ScaleFactor)
    ScaleFactor.radius = 1;
    ScaleFactor.velocity = 1;
end
if nargin < 4 || isempty(Perturbation)
    Perturbation.density = 0;
    Perturbation.C_D = 0;
    Perturbation.C_L = 0;
end
%Vehicle Parameters, these could be passed in instead
S = 15.8; %reference wing area, m^2
m = 2804; %mass, kg

%State Variables:
% x = [theta phi r psi gamma V]

s = 1;
theta = x(s+1);   %longitude, rad
phi = x(s+2);     %latitude, rad
r = x(s);       %vehicle radius, m
psi = x(s+5);     %heading angle, rad, 0 -> due East
gamma = x(s+4);   %flight path angle, rad
V = x(s+3);       %velocity, m/s

g = mu/r^2/ScaleFactor.velocity^2/ScaleFactor.radius;      % gravity, m/s^2
h = r*ScaleFactor.radius-r_eq;                             % Altitude
[rho,a,hs] = MarsAtmosphericDensity(h);                       % Density and Speed of Sound
rho = rho+Perturbation.density;
M = V*ScaleFactor.velocity/a;                              % Mach Number
[C_D,C_L] = AerodynamicCoefficients(M); 
C_D = C_D + Perturbation.C_D;
C_L = C_L + Perturbation.C_L;
temp = .5*rho*V^2*S/m*ScaleFactor.radius;
L = temp*C_L;
D = temp*C_D;

r_dot = V*sin(gamma);
V_dot = -D-g*sin(gamma);
theta_dot = V/r*cos(gamma)*cos(psi)/cos(phi);
gamma_dot = (L*cos(sigma) - (g-V^2/r)*cos(gamma))/V;% + 2*omega_p*cos(psi)*cos(phi);
phi_dot = V/r*cos(gamma)*sin(psi);
psi_dot = 1/(V*cos(gamma))*(L*sin(sigma)-V^2/r*cos(gamma)^2*cos(psi)*tan(phi));% + 2*omega_p*(tan(gamma)*sin(psi)*cos(phi)-sin(phi));

dX = [r_dot; theta_dot; phi_dot;  V_dot; gamma_dot; psi_dot ];
end
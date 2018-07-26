% FINALLATLON Computes the coordinates of a point with given downrange and
% crossrange from an initial point.
%   FINALLATLON(THETA, PHI, PSI, DR, CR) computes the desired latitude and
%   longitude given an initial longitude THETA, latitude PHI, heading angle
%   PSI, and downrange and crossrange requirements, DR and CR respectively.
%
%   The angle inputs must be in radians and the distance inputs must be in
%   kilometers.
%
%   The radius of Mars is hard-coded, so don't use it for Earth!
function [theta_f,phi_f] = FinalLatLon(theta,phi,psi,DR,CR)

r_p = 3397; %km
LF = acos(cos(DR/r_p)*cos(CR/r_p));
zeta = asin(sin(CR/r_p)/sin(LF));
phi_f = -asin(cos(zeta-psi+pi/2)*cos(phi)*sin(LF)+sin(phi)*cos(LF));
theta_f = theta + asin(sin(zeta-psi+pi/2)*sin(LF)/cos(phi_f));

end


% RANGE Computes the downrange and crossrange of an entry vehicle given the
% initial lat/long and heading, and the current lat/long. All inputs should
% be in radians. The outputs are given in kilometers. The function is 
% currently hard coded for Mars; in the future, the planet radius ( or
% planet structure) should be a required input as well.
%
% Note:
% The output is slightly buggy when the current and initial coordinates are 
% equal or nearly equal. That is why the real/isnan fixes are in place. The
% equations are undefined at this point due to a divide by zero error.

function [DR,CR] = Range(theta0,phi0,psi0,theta,phi)

r_p = 3397; %km
LF = real(acos(sin(phi)*sin(phi0)+cos(phi).*cos(phi0).*cos(theta-theta0)));
sig = asin(sin(theta-theta0).*cos(phi)./sin(LF));
zeta = sig+psi0-pi/2;

DR = real(r_p*acos(cos(LF)./cos(asin(sin(LF).*sin(zeta)))));
CR = -real(r_p*asin(sin(LF).*sin(zeta)));
DR(isnan(DR)) = 0;
CR(isnan(CR)) = 0;

end
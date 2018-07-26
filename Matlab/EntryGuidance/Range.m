% RANGE Computes the downrange and crossrange of an entry vehicle given the
% initial lat/long and heading, and the current lat/long. All inputs should
% be in radians. The outputs are given in kilometers. The function is
% currently hard coded for Mars; in the future, the planet radius ( or
% planet structure) should be a required input as well.

function [DR,CR] = Range(lon0,lat0,heading0,lonc,latc)

r_p = 3397; %km

d13 = acos(sin(latc).*sin(lat0)+cos(latc).*cos(lat0).*cos(lonc-lon0));
% if abs(d13) < 1e-4
%     DR=0.*d13;
%     CR=0.*d13;
% else
    psi12 = heading0;
    PHI = sign(lonc-lon0).*acos( (sin(latc) - sin(lat0).*cos(d13))./(cos(lat0).*sin(d13)) );
    psi13 = pi/2 - PHI;
    CR = asin(sin(d13).*sin(psi12-psi13));
    DR = r_p*acos(cos(d13)./cos(CR));
    CR = -r_p*CR;
    
% end
end
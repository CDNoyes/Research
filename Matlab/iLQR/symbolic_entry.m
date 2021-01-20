clear
clc
syms h rp v fpa mu rho0 hs cd cl S m u s hscale fscale sscale r D L g rho

scale = [hscale, fscale, sscale];
x = [h; fpa; s];
xu = [x;u];

% r = h + rp;
% g = mu/r^2;
% rho = rho0*exp(-h/hs);
% f = 0.5*rho*v^2 * S/m;
% D = f*cd;
% L = f*cl;


hdot = v*sin(fpa); 
sdot = v*cos(fpa);
vdot = -D - g*sin(fpa);
fdot = L/v*u + (v/r-g/v)*cos(fpa);

dvdx = jacobian(vdot, [h, fpa, s, u]);


hprime = hdot/vdot;
Jh = jacobian(hdot, [h, fpa, s, u]);
Jhp = jacobian(hprime, [h, fpa, s, u]);
Jhp2 = Jh/vdot - hdot*dvdx/(vdot^2);

xdot = [vdot, hdot, fdot, sdot].'; % zero is for the control term 

J = jacobian(xdot, [h, fpa, s, u]); % jacobian in time
% Jv = jacobian(vdot, v);
Jp = jacobian(xdot/vdot, xu); % jacobian in velocity
% Jp2 =  J/vdot - (repmat(xdot,1,4).*repmat(dvdx, 3,1))/(vdot^2); %
% equivalent

function [rho, a, hs,a_prime,a_dot,rho_dot,rho_ddot,rho_dddot] = MarsAtmosphericDensity(h,h_dot)
% input  h, in meters
% input  h_dot, in meters/second (recall h_dot == r_dot)
% output rho, in kg/m^3
% output a, local speed of sound, m/s
% output hs, scale height, m
% output a_dot = m/s^2

if nargin == 1
    h_dot = 0;
end

%Exponential Curve Fit
% Density Calculation:
% beta = [-4.324,-9.204e-5, -1.936e-11, -7.507e-15, 4.195e-20];
% rho = exp(beta(5)*h^4 + beta(4)*h^3 + beta(3)*h^2 + beta(2)*h + beta(1));
% % "Scale Height"
% hs = -1/(beta(2)+2*beta(3)*h+3*beta(4)*h^2+4*beta(5)*h^3);

hs = 9354.5;
rho = 0.0158*exp(-h/hs); %Simple exponential model
rho_dot = rho*h_dot/-hs;
rho_ddot = rho*(h_dot/-hs)^2;
rho_dddot = rho*(h_dot/-hs)^3;

% Local Speed of Sound
b = [223.8 -0.2004e-3 -1.588e-8 1.404e-13];
a = b(4)*h.^3 + b(3)*h.^2 + b(2)*h + b(1);
a_prime = (3*b(4)*h.^2 + 2*b(3).*h + b(2));
a_dot = a_prime.*h_dot;

end
function dx = Coriolis(x,planetModel)

omega_p = planetModel.omega;
phi = x(3);     %latitude, rad
gamma = x(5);   %flight path angle, rad
psi = x(6);     %heading angle, rad, 0 -> due East

% r_dot = 0;
% V_dot = 0;
% theta_dot = 0;
gamma_dot = 2*omega_p*cos(psi)*cos(phi);
% phi_dot = 0;
psi_dot = 2*omega_p*(tan(gamma)*sin(psi)*cos(phi)-sin(phi));

dx = [0;0;0;gamma_dot;0;psi_dot];

end
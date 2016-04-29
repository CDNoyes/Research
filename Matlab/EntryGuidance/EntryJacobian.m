function [J,B] = EntryJacobian(x,sigma,planetModel,vehicleModel,ScaleFactor)
if nargin < 5 || isempty(ScaleFactor)
    ScaleFactor.radius = 1;
    ScaleFactor.velocity = 1;
end

omega_p = 0; %planetModel.omega;

% theta = x(2);
phi = x(:,3);
r = x(:,1);
gamma = x(:,5);
psi = x(:,6);
v = x(:,4);

[g,L,D,hs,M,a,rho,rho_dot,q,C_D,C_D_dot,C_L,dadr,C_D_prime,C_L_prime] = EntryForces(x,planetModel,vehicleModel,ScaleFactor);

dMda = -M/(a/ScaleFactor.velocity); 
dMdr = dMda.*dadr*ScaleFactor.velocity; % dadr is our model of the speed of sound as a function of altitude (because temperature varies with altitude)
dCLdr = C_L_prime.*dMdr;
dCDdr = C_D_prime.*dMdr; % Prime is wrt Mach for CD and CL
dCLdv = C_L_prime./(a/ScaleFactor.velocity);
dCDdv =  C_D_prime./(a/ScaleFactor.velocity);

J = zeros(6,6,length(r));

%Partials of theta_dot wrt
J(2,3,:) = v.*cos(gamma).*cos(psi).*sin(phi)./r./cos(phi).^2;   % phi
J(2,1,:) = -v.*cos(gamma).*cos(psi)./r.^2./cos(phi);           % r
J(2,5,:) = -v.*sin(gamma).*cos(psi)./r./cos(phi);             % gamma
J(2,6,:) = -v.*cos(gamma).*sin(psi)./r./cos(phi);              % psi
J(2,4,:) = cos(gamma).*cos(psi)./r./cos(phi);                % v

%Partials of phi_dot wrt
J(3,1,:) = -v.*cos(gamma).*sin(psi)./r.^2;    % r
J(3,5,:) = -v.*sin(gamma).*sin(psi)./r;      % gamma
J(3,6,:) = v.*cos(gamma).*cos(psi)./r;      % psi
J(3,4,:) = cos(gamma).*sin(psi)./r;         % v

%Partials of r_dot wrt 
J(1,5,:) = v.*cos(gamma);  % gamma
J(1,4,:) = sin(gamma);    % V

%Partials of gamma_dot wrt
J(5,3,:) = -2*omega_p*cos(psi).*sin(phi);  % phi
J(5,1,:) = -v.*cos(gamma)./r.^2+(-L./hs).*cos(sigma)./v + 2.*g.*cos(gamma)./r./v + cos(sigma).*L.*dCLdr./C_L./v; % r
J(5,5,:) = (g./v-v./r).*sin(gamma);          % gamma
J(5,6,:) = -2.*omega_p.*sin(psi).*cos(phi); % psi
J(5,4,:) = L.*cos(sigma)./v.^2 + (g./v.^2+1./r).*cos(gamma) + cos(sigma).*q.*dCLdv./v; %v

%Partials of psi_dot wrt
J(6,3,:) = -v.*cos(psi).*cos(gamma)./r./cos(phi).^2 + 2.*omega_p.*(-tan(gamma).*sin(psi).*sin(phi)-cos(phi)); % phi
J(6,1,:) = (L./hs).*sin(sigma)./v./cos(gamma) + v.*cos(psi).*cos(gamma).*tan(phi)./r.^2 - sin(sigma).*q.*dCLdr./v./cos(gamma); % r
J(6,5,:) = v.*sin(gamma).*cos(psi).*tan(phi)./r + L.*sin(sigma).*sin(gamma)./v./cos(gamma).^2 + 2.*omega_p.*(sin(psi).*cos(phi)./cos(gamma).^2); % gamma
J(6,6,:) = v.*sin(psi).*cos(gamma).*tan(phi)./r + 2*omega_p*tan(gamma).*cos(psi).*cos(phi); % psi
J(6,4,:) = -cos(psi).*cos(gamma).*tan(phi)./r + L.*sin(sigma)./v.^2./cos(gamma) - sin(sigma).*L.*dCLdv./C_L./v./cos(gamma); % v

%Partials of v_dot wrt
J(4,1,:) = 2*g.*sin(gamma)./r + D./hs - q.*dCDdr;          % r
J(4,4,:) = -2*D./v - q.*dCDdv;                             % v
J(4,5,:) = -g.*cos(gamma);                                 % gamma

%Partials wrt to sigma
B = zeros(6,1,length(r));
B(5,1,:) = -L./v.*sin(sigma);
B(6,1,:) = -L./v./cos(gamma).*cos(sigma);

end
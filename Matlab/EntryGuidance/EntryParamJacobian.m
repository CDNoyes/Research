function J = EntryParamJacobian(x,sigma,planetModel,vehicleModel,ScaleFactor)
if nargin < 5 || isempty(ScaleFactor)
    ScaleFactor.radius = 1;
    ScaleFactor.velocity = 1;
end

omega_p = 0; %planetModel.omega;

% theta = x(2);
% phi = x(:,3);
r = x(:,1);
gamma = x(:,5);
% psi = x(:,6);
v = x(:,4);

[g,L,D,hs,M,a,rho,rho_dot,q,C_D,C_D_dot,C_L] = EntryForces(x,planetModel,vehicleModel,ScaleFactor);

% dMda = -M/(a/ScaleFactor.velocity); 
% dMdr = dMda.*dadr*ScaleFactor.velocity; % dadr is our model of the speed of sound as a function of altitude (because temperature varies with altitude)
% dCLdr = C_L_prime.*dMdr;
% dCDdr = C_D_prime.*dMdr; % Prime is wrt Mach for CD and CL
% dCLdv = C_L_prime./(a/ScaleFactor.velocity);
% dCDdv =  C_D_prime./(a/ScaleFactor.velocity);

J = zeros(6,3,length(r));

%Partials of theta_dot wrt
%All zeros


%Partials of phi_dot wrt
% All zeros

%Partials of r_dot wrt 
% All zeros


%Partials of gamma_dot wrt
J(5,2,:) = -sin(sigma)./v./cos(gamma).*L./C_L;         % CL
J(5,3,:) = -sin(sigma)./v./cos(gamma).*L./rho;   % rho


%Partials of psi_dot wrt
J(6,2,:) = cos(sigma)./v.*L./C_L; % CL
J(6,3,:) = cos(sigma)./v.*L./rho; % rho


%Partials of v_dot wrt
J(4,1,:) = -D./C_D;          % CD
% J(4,2,:) = 0;              % CL
J(4,3,:) = -D./rho;          % rho

end
function DragSymbolic()

syms S m CD V rho rho_dot rho_ddot V_dot V_ddot CD_dot g gamma gamma_dot u r D(rho,V,CD) L(rho,V)

x = [rho; V; CD]; % The time-varying elements needed for the jacobian
dxdt = [rho_dot; V_dot; CD_dot];
c = 0.5*S/m; % Constant stuff, left off for easier readability

% First Derivative:
Ddot = jacobian(D,x)*dxdt;
Ddot = subs(Ddot,[diff(D(rho, V, CD), rho),diff(D(rho, V, CD), V),diff(D(rho, V, CD), CD)], [D/rho,2*D/V,D/CD]);
%Check my handwritten formulation
D_dot = D*(rho_dot/rho + 2*V_dot/V+CD_dot/CD);
assert(vpa(simple(Ddot-D_dot))==0)

% Second Derivative:                           
X = [x;dxdt];
CD_ddot = 0; %Add to syms to remove assumption
dXdt = [dxdt; rho_ddot; V_ddot; CD_ddot]; %CD_ddot = 0 assumed here
Dddot = jacobian(D_dot,X)*dXdt;
Dddot = subs(simple(Dddot),[diff(D(rho, V, CD), rho),diff(D(rho, V, CD), V),diff(D(rho, V, CD), CD),CD*V^2*rho], [D/rho,2*D/V,D/CD,D]);
%Check my formula again
D_ddot = D_dot^2/D + D*(rho_ddot/rho-(rho_dot/rho)^2 + 2*(V_ddot/V-(V_dot/V)^2) + (CD_ddot/CD-(CD_dot/CD)^2));
simple(Dddot-D_ddot)
assert(vpa(simple(Dddot-D_ddot))==0)

%% D_ddot = a + bu
% gamma_dot = L*u/V + (V/r-g/V)*cos(gamma);
% g_dot = -2*g*V*sin(gamma)/r;
% V_ddot = -Ddot-g_dot*sin(gamma)-g*cos(gamma)*gamma_dot;


end
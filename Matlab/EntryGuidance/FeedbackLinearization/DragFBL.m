% DRAGFBL Returns the parameters a,b such that the second time rate of
% change of drag = a+bu.
%   Even when using a control law not based on FBL, this function is used
%   in the sliding mode observer

function [a,b,ahat,bhat] = DragFBL(g,L,D,r,V,gamma,rho,rho_dot,D_dot)

V_dot = -D-g*sin(gamma);
g_dot = -2*g*V*sin(gamma)/r;
h_dot = V*sin(gamma);


if isempty(D_dot) %When using an observer, we used the observed estimate instead of the model estimate
    D_dot = D*(rho_dot/rho + 2*V_dot/V); %+ D*CD_dot/CD %neglect the variation in C_D
end

E_dot = -V*D;
E_ddot = -V_Dot*D - V*D_dot;

a1 = D_dot*(rho_dot/rho + 2*V_dot/V) - D*(rho_dot^2/rho^2 + 2*V_dot^2/V^2);
a2 = -2*D/V*(D_dot+g_dot*sin(gamma));
a3 = -2*D*g*cos(gamma)^2 * (1/r - g/V^2);
a4 = -D*(rho_dot/rho)^2 + D*rho_dot/rho/h_dot*(-g-D*sin(gamma)+V^2/r*cos(gamma)^2);
a = a1 + a2 + a3 + a4;

b1 = -2*D*L*g*cos(gamma)/V^2;
b2 = D*L/h_dot*rho_dot/rho*cos(gamma);
b = b1+b2;

% The parameters with respect to Energy instead of Time
ahat = a/(E_dot^2) - D_dot*E_ddot/(E_dot^3);
bhat = b/(E_dot^2);

end
function [a,b] = DragFBL(g,L,D,r,V,gamma,rho,rho_dot,D_dot)

if isempty(D_dot)
    D_dot = -D*V*sin(gamma)/hs - 2*D/V*(D+g*sin(gamma));
end
% 
% a = -D_dot*V*sin(gamma)/hs + D*(D+g*sin(gamma))*sin(gamma)/hs + ...
%     D*cos(gamma)^2*(g-V^2/r)/hs - 2*D_dot*(D+g*sin(gamma))/V^2 + ... %the second term has a minus
%     -2*D/V^2*(D+g*sin(gamma))^2 - 2*D*D_dot/V + 4*D*g*sin(gamma)^2/r + ...
%     2*D/V^2*g*cos(gamma)^2*(g-V^2/r);
% 
% b = -D*L*cos(gamma)*(2*g/V^2 +(1/hs));


V_dot = -D-g*sin(gamma);
g_dot = -2*g*V*sin(gamma)/r;
h_dot = V*sin(gamma);

a1 = D_dot*(rho_dot/rho + 2*V_dot/V) - D*(rho_dot^2/rho^2 + 2*V_dot^2/V^2);
a2 = -2*D/V*(D_dot+g_dot*sin(gamma));
a3 = -2*D*g*cos(gamma)^2 * (1/r - g/V^2);
a4 = -D*(rho_dot/rho)^2 + D*rho_dot/rho/h_dot*(-g-D*sin(gamma)+V^2/r*cos(gamma)^2);
a = a1 + a2 + a3 + a4;

b1 = -2*D*L*g*cos(gamma)/V^2;
b2 = D*L/h_dot*rho_dot/rho*cos(gamma);
b = b1+b2;

end
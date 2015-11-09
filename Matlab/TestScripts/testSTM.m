% test the STM methods


%% Vanderpol oscillator tests

f = @(t,x,u,p) [x(2);-x(1)+p*(1-x(1).^2).*x(2) + u(t) ];
jac = @(x,u,p) [0, 1; -1-2*p*x(1).*x(2)  p*(1-x(1).^2)];

u = @(t) 0;
% Try out the ControlWrapper functionality

%% Vanderpol oscillator tests
clc; clear;
f = @(t,x,u,p) [x(2);-x(1)+p*(1-x(1).^2).*x(2) + ControlWrapper(t,u)];
jac = @(x,u,p) [0, 1; -1-2*p*x(1).*x(2)  p*(1-x(1).^2)];
dim.state = 2;
dim.control = 1;
dim.param = 1;
index.state = 1:dim.state;
index.control = dim.state+1:dim.state+dim.control;
index.param = index.control(end)+1:index.control(end)+dim.param;
u = @(t) cos(t)+sin(t);

%By using the wrapped control, I can both take the Jacobian numerically and
%integrate the equations using a control law

x0 = [10;2;zeros(dim.control,1);0.25];
J0 = ComplexDiff(@(X) f(0,X(index.state),X(index.control),X(index.param)),x0);
J0_ana = jac(x0(index.state),x0(index.control),x0(index.param));
[t,x] = ode45(@(t,x) f(t,x,u,0.25),[0,12],x0(index.state));
plot(x(:,1),x(:,2))
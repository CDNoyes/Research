% Test the numerical differentiation utilities

%% Scalar function of 2 variables with constant Hessian
clear; clc;
H = [1, 2;2,-0.5];
testFunction = @(x) 0.5*(x.'*H*x);
exactGradient = @(x) x.'*H;
exactHessian = @(x) H;

testPoint = [1;1];
complexError = norm(exactGradient(testPoint) - ComplexDiff(testFunction,testPoint));
centralError = norm(exactGradient(testPoint) - CentralDiff(testFunction,testPoint));
forwardError = norm(exactGradient(testPoint) - ForwardDiff(testFunction,testPoint));
hessianError = norm(exactHessian(testPoint)  - Hessian(testFunction,testPoint));

%% Scalar function of 2 variables, non-constant Hessian:
clear; clc;
testFunction = @(x) x(1)*exp(x(2));
exactGradient = @(x) [exp(x(2)), x(1)*exp(x(2))];
exactHessian = @(x) [0, exp(x(2)); exp(x(2)), x(1)*exp(x(2))];

testPoint = [2;3];

complexError = norm(exactGradient(testPoint) - ComplexDiff(testFunction,testPoint));
centralError = norm(exactGradient(testPoint) - CentralDiff(testFunction,testPoint));
forwardError = norm(exactGradient(testPoint) - ForwardDiff(testFunction,testPoint));
hessianError = norm(exactHessian(testPoint)  - Hessian(testFunction,testPoint));

%% Vector-valued function of 2 variables:
clear; clc;
testFunction = @(x) [x(2);-x(1)+0.25*(1-x(1)^2)*x(2)];
exactJacobian = @(x) [0, 1; -1-2*0.25*x(1)*x(2)  0.25*(1-x(1)^2) ];

testPoint = [1;-2];
x = testPoint;
exactHessian(1:2,1:2,1) = zeros(2);
exactHessian(1:2,1:2,2) = [-0.5*x(2),-0.5*x(1); -0.5*x(1),0];

complexError = norm(exactJacobian(testPoint) - ComplexDiff(testFunction,testPoint));
centralError = norm(exactJacobian(testPoint) - CentralDiff(testFunction,testPoint));
forwardError = norm(exactJacobian(testPoint) - ForwardDiff(testFunction,testPoint));
numHessian = Hessian(testFunction,testPoint);
hessianError = norm([norm(exactHessian(:,:,1)-numHessian(:,:,1)),norm(exactHessian(:,:,2)-numHessian(:,:,2))]);

%% VDP Oscillator: vector valued with two states, one control, and one parameter
clear; clc;
[F,J,H] = VDP();
[dim,ind] = Dimension(2,1,1);

x0 = [3;5];
u0 = 0;
mu = 0.25;
X0 = [x0;u0;mu];

%Analytical Results
f = F(0,x0,u0,mu);
j = J(x0,u0,mu);
h = H(x0,u0,mu);

%Numerical Results
F1 = @(X) F(0,X(ind.state),X(ind.control),X(ind.parameter));
jn = ComplexDiff(F1,X0);
hn = Hessian(F1,X0);
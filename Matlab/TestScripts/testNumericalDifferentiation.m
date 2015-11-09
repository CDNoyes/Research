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
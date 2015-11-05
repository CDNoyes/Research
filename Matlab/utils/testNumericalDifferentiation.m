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

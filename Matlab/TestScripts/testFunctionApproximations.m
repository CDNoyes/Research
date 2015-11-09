% test the linearization and quadratization functions

%% Vanderpol oscillator
clear; clc;

testFunction = @(x) [x(2);-x(1)+0.25*(1-x(1).^2).*x(2)];
exactJacobian = @(x) [0, 1; -1-2*0.25*x(1).*x(2)  0.25*(1-x(1).^2) ];

appPoint = [1;-2]; %Point of approximation
fLinear = Linearize(testFunction,appPoint,[]);%exactJacobian);
eLinear = @(x) norm(fLinear(x)-testFunction(x)); %Error between the linearized function and the true function

fQuad = Quadratize(testFunction,appPoint,exactJacobian);
eQuad = @(x) norm(fQuad(x)-testFunction(x));
for k = 1%:500

[x1,x2] = ndgrid(0:.05:2,-3:.05:-1);
clear eL eQ
tLin = 0;
tQuad = 0;
for i = 1:size(x1,1)
    for j = 1:size(x1,2)
        tic;
        eL(i,j) = eLinear([x1(i,j);x2(i,j)]);
        tLin = tLin + toc;
        tic;
        eQ(i,j) = eQuad([x1(i,j);x2(i,j)]);
        tQuad = tQuad + toc;
    end
end

ratio(k) = tQuad/tLin;
end
figure
surf(x1,x2,eL)
xlabel('x_1')
ylabel('x_2')
zlabel('||error||')

figure
surf(x1,x2,eQ)
xlabel('x_1')
ylabel('x_2')
zlabel('||error||')
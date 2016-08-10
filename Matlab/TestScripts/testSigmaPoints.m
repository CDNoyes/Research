% Test SigmaPoints

clear;clc;

n = 2;
x = zeros(n,1);
% R = Regularize(eye(n) + -.1 + 0.2*rand(n));
R = [1,0.9;0.9,1];
% R = diag([5,1]);
[X,W] = SigmaPoints(x,R);


[xMean,P] = SigmaEval(X,W);

errState = abs(x-xMean);
errCov = abs(P-R); 

xn = mvnrnd(x,R,1000);
pn = mvnpdf(xn,x',R);

nDev = [1,2,3];


figure
% plot(x(1),x(2),'rx','MarkerSize',10)
scatter(xn(:,1),xn(:,2),8,pn)
hold all
plot(X(1,:),X(2,:),'k*','MarkerSize',10)
DrawNormalEllipse(x,R,nDev);
% plot(ell(1,:),ell(2,:),'k--','LineWidth',4)

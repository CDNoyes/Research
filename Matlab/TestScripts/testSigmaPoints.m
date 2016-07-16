% Test SigmaPoints

clear;

n = 3;
x = zeros(n,1);
R = Regularize(eye(n) + -.1 + 0.2*rand(n));
[X,W] = SigmaPoints(x,R);


xMean = sum(repmat(W.state,n,1).*X,2);
Xd = X-repmat(xMean,1,size(X,2));
P = (repmat(W.cov,n,1).*Xd)*Xd';

errState = abs(x-xMean);
errCov = abs(P-R); 
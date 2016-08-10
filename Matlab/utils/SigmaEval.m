% SIGMAEVAL(X,W) Estimates the mean and covariance of an unscented
% transform.
%   SIGMAEVAL(X,W) Uses the sigma points in X (after having undergone some
%   nonlinear transformation) and the weights W to estimate the 

function [xMean,P] = SigmaEval(X,W)

[n,m] = size(X);

xMean = sum(repmat(W.state,n,1).*X,2);
Xd = X-repmat(xMean,1,m);
P = (repmat(W.cov,n,1).*Xd)*Xd';



end
% SIGMAPOINTS Computes the sigma points of an unscented transform.
%   [X,W] = SIGMAPOINTS(XMEAN,XCOV) Given a mean point XMEAN, and a
%   corresponding covariance matrix XCOV, computes the (2n+1) sigma points
%   stored columnwise in output X for use in an unscented transform. The
%   weights W are used to evaluate the mean and covariance using the sigma
%   points after a nonlinear transformation.
%
%   The values of alpha, beta, and k used here are not the only choice, but
%   have been selected based on the recommendations of the following reference:
%   http://nbviewer.jupyter.org/github/sbitzer/UKF-exposed/blob/master/UKF.ipynb

function [X,W] = SigmaPoints(xMean,xCov)

L = length(xMean);
xMean = xMean(:);

alpha = 1;
beta = 0; % This value depends on the distribution of x, for gaussians 2 is optimal
k = L;
lambda = alpha^2*(L+k)-L;
W.state = [lambda/(L+lambda), ones(1,2*L)/2/(L+lambda)];
W.cov = [lambda/(L+lambda) + (1-alpha^2+beta), ones(1,2*L)/2/(L+lambda)];
R = chol((L+lambda)*xCov,'lower');
X = [xMean, repmat(xMean,1,L)+R, repmat(xMean,1,L)-R]; %Each column of X is a sigma pt

end
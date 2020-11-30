function [X,W] = UnscentedTransform(mean, covariance, scale)
% Computes sigma points X and weights W 

if nargin < 3
   scale = 3; 
end

n = max(size(mean));

S = chol((n+scale)*covariance);
X = [mean; pp(mean, S); pp(mean,-S)]';
W = ones(1,2*n+1)*0.5/(n+scale);
W(1) = scale/(n+scale);


function c = pp(a,b)
c = bsxfun(@plus,a,b);
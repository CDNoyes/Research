function DrawUniformBounds(xMean,xCov,lineSpec)

if nargin < 3
   lineSpec = '--';
end

A = chol(3*Regularize(xCov,1e-12),'lower');
X = repmat(xMean(:),1,4);
% y = eye(2);
y = [1 -1;1 1];
% ell = X + [sqrt(3)*A*y,-sqrt(3)*A*y];
ell = X + [A*y,-A*y];

plot(ell(1,1:2),ell(2,1:2),lineSpec)
plot(ell(1,2:3),ell(2,2:3),lineSpec)
plot(ell(1,3:4),ell(2,3:4),lineSpec)
plot(ell(1,[4,1]),ell(2,[4,1]),lineSpec)



end
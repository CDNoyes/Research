function DrawNormalEllipse(xMean,xCov,nDev,lineSpec)

if nargin < 3
   nDev = 3;
   lineSpec = '--';
end
if nargin < 4
   lineSpec = '--'; 
end

theta = linspace(0,2*pi);
xy = [cos(theta); sin(theta)];
A = chol(Regularize(xCov,1e-12),'lower');
X = repmat(xMean(:),1,100);

for i = 1:length(nDev)
    ell = X + nDev(i)*A*xy;
    plot(ell(1,:),ell(2,:),lineSpec)
end
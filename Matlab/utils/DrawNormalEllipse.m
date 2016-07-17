function DrawNormalEllipse(xMean,xCov,nDev)

if nargin < 3
   nDev = 3;
end

theta = linspace(0,2*pi);
xy = [cos(theta)+xMean(1); sin(theta)+xMean(2)];
A = chol(xCov,'lower');
ell = A*xy;

for i = 1:length(nDev)
    plot(nDev(i)*ell(1,:),nDev(i)*ell(2,:))
end
%OPTIMIZEBANKANGLEDYNAMICS Optimizes the gains and prediction horizon used
%in the second-order bank angle dynamics.
%   OPTIMIZEBANKANGLEDYNAMICS(GAINS,SIGMA) computes the error between a
%   bank angle profile given by the function handle SIGMA and a bank angle
%   profile resulting from a constrained second-order system using the
%   parameters in GAINS.
function err = optimizeBankAngleDynamics(gains,Sigma)
tf = 265.0143;
dtr = pi/180;
sigmaMin = 18.19*dtr;
sigmaMax = 87.13*dtr;
rateMax = 20*dtr;
accMax = 5*dtr;
lim.rate = rateMax;
lim.acceleration = accMax;
lim.angleMax = sigmaMax;
lim.angleMin = sigmaMin;
[T,S] = ode45(@BankAngleDynamics,[0,tf],[Sigma(0);0],[],Sigma,lim,gains);
err = norm(Sigma(T)'/dtr-Saturate(S(:,1),-sigmaMax,sigmaMax)/dtr);

end
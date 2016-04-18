%OPTIMIZEBANKANGLEDYNAMICS Optimizes the gains and prediction horizon used
%in the second-order bank angle dynamics.
%   OPTIMIZEBANKANGLEDYNAMICS(GAINS,SIGMA,LIM) computes the error between a
%   bank angle profile given by the function handle SIGMA and a bank angle
%   profile resulting from a second-order system constrained by the limits
%   set in LIM using the parameters in GAINS.
function err = optimizeBankAngleDynamics(gains,Sigma,lim,tf)

dtr = pi/180;
gains(3) = 0;
[T,S] = ode45(@BankAngleDynamics,[0,tf],[Sigma(0);0],[],Sigma,lim,gains);
% err = norm(Sigma(T)'/dtr-Saturate(S(:,1),-lim.angleMax,lim.angleMax)/dtr);
err = trapz(T,(Sigma(T)'/dtr-Saturate(S(:,1),-lim.angleMax,lim.angleMax)/dtr).^2);
% err = sum((Sigma(T)'/dtr-Saturate(S(:,1),-lim.angleMax,lim.angleMax)/dtr).^2);

end
%LPF Implements a low pass filter in the time domain.
%   LPF(TIME,SIGNAL,TIMECONSTANT) Computes the low-pass filtered value of
%   SIGNAL at the time(s) in the vector TIME with frequency (in Hz) equal
%   to TIMECONSTANT^-1. SIGNAL should be a function handle taking only a
%   single argument.

function LPsignal = LPF(time,signal,timeConstant)

f = @(t,x) (-x+signal(t))/timeConstant;
% opt = odeset('AbsTol',1e-9,'RelTol',1e-9);
opt = odeset('AbsTol',1e-12,'RelTol',1e-12);
if isscalar(time)
    [~,X] = ode45(f,[0,time],signal(0),opt);
    LPsignal = X(end);
else
    [~,LPsignal] = ode45(f,time,signal(time(1)),opt);

end



end
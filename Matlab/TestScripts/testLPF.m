clear;
t = linspace(0,10,1000);
tn = linspace(t(1),t(end),100);
n = -.3+.6*rand(1,length(tn));
s = @(T) T+sin(length(tn)*T)+interp1(tn,n,T,'spline');
tau = .1;

sLP = LPF(t,s,tau);

plot(t,s(t),'r',t,sLP,'k')
legend(['Noisy Signal(',num2str(1/(t(end)/length(tn))),' Hz)'],['Filtered Signal (',num2str(1/tau),' Hz)'])
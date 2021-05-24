function [hm, fm, sm, hs, fs, ss]=get_stats(sol)
hm = sol.mean(1,end)/1000;
fm = sol.mean(2,end)*180/pi;
sm = sol.mean(3,end);
if sm(end) > 100e3
    k = 1000;
else
    k = 1;
end
sm = sm/k;

hv = sol.var(1,end)/1000^2;
fv = sol.var(2,end)*(180/pi)^2;
sv = sol.var(3,end)/k^2;

hs = hv^0.5;
fs = fv.^0.5;
ss = sv.^0.5;

disp(['hf  = ',num2str(hm),'  +/- 3*',num2str(hv.^0.5),' km (3\sigma low = ',num2str(hm-3*hs),')'])
disp(['fpa = ',num2str(fm),'  +/- 3*',num2str(fv.^0.5),' deg'])
disp(['sf  = ',num2str(sm),'  +/- 3*', num2str(sv.^0.5),' km (3\sigma = ',num2str(3*ss),')'])
disp(' ')
function EntryPlots(t,x)
dtr = pi/180;
r_eq = 3397e3;
hkm = (x(:,1)-r_eq)/1000;
% Altitude vs Velocity
plot(x(:,4),hkm)
xlabel('Velocity (m/s)')
ylabel('Altitude (km)')
title(['Final altitude: ',num2str(hkm(end)), ' km'])
%DR v CR
[DR,CR] = Range(x(1,2),x(1,3),x(1,6),x(:,2),x(:,3));

figure
plot(CR,DR)
ylabel('Downrange (km)')
xlabel('Crossrange (km)')
title(['Downrange: ',num2str(DR(end)),' km, Crossrange: ',num2str(CR(end))])

figure
plot(x(:,3)/dtr,x(:,2)/dtr)
xlabel('Latitude (deg)')
ylabel('Longitude (deg)')
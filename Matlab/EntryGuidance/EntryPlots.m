function EntryPlots(t,x)
dtr = pi/180;
planet = Mars();
r_eq = planet.radiusEquatorial;
hkm = (x(:,1)-r_eq)/1000;

% Altitude vs Velocity
figure
ParachuteDeploymentConstraints(true);
hold on
plot(x(:,4),hkm)
xlabel('Velocity (m/s)')
ylabel('Altitude (km)')
title(['Final altitude: ',num2str(hkm(end)), ' km'])
disp(['Final altitude: ',num2str(hkm(end)), ' km'])

%DR v CR
[DR,CR] = Range(x(1,2),x(1,3),x(1,6),x(:,2),x(:,3));
figure
plot(CR,DR)
ylabel('Downrange (km)')
xlabel('Crossrange (km)')
title(['Downrange: ',num2str(DR(end)),' km, Crossrange: ',num2str(CR(end)), ' km'])

%Lat/Lon
figure
plot(x(1,3)/dtr,x(1,2)/dtr,'ko')
hold on
plot(x(:,3)/dtr,x(:,2)/dtr)
xlabel('Latitude (deg)')
ylabel('Longitude (deg)')
legend('Initial point','Location','Best')

%Flight path angle
figure
plot(t,x(:,5)/dtr)
xlabel('Time (s)')
ylabel('FPA (deg)')
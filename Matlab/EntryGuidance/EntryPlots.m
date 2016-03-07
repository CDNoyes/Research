%ENTRYPLOTS Plot some standard graphs useful in examining entry
%trajectories.

function EntryPlots(trajectorySummary)
dtr = pi/180;
planet = Mars();
ts = trajectorySummary;
t = ts.time;
x = ts.state;
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

%DR v CR
figure
plot(ts.CR,ts.DR)
ylabel('Downrange (km)')
xlabel('Crossrange (km)')
title(['Downrange: ',num2str(ts.DR(end)),' km, Crossrange: ',num2str(ts.CR(end)), ' km'])

%Lat/Lon
figure
plot(x(1,3)/dtr,x(1,2)/dtr,'ko')
hold on
plot(ts.target.lat/dtr,ts.target.lon/dtr,'kx')
plot(x(:,3)/dtr,x(:,2)/dtr)
xlabel('Latitude (deg)')
ylabel('Longitude (deg)')
legend('Initial point','Target','Location','Best')

%Flight path angle
figure
subplot 211
plot(t,x(:,5)/dtr)
xlabel('Time (s)')
ylabel('FPA (deg)')
subplot 212
plot(ts.energy_norm,x(:,5)/dtr)
xlabel('Normalized Energy')
ylabel('FPA (deg)')

%Heading angle
figure
subplot 211
plot(t,x(:,6)/dtr)
xlabel('Time (s)')
ylabel('Heading (deg)')
subplot 212
plot(ts.energy_norm,x(:,6)/dtr)
xlabel('Normalized Energy')
ylabel('Heading (deg)')
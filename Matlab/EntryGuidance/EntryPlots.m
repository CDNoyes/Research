%ENTRYPLOTS Plot some standard graphs useful in examining entry
%trajectories.

function output = EntryPlots(t,x,DR_t,CR_t)
dtr = pi/180;
planet = Mars();
r_eq = planet.radiusEquatorial;
hkm = (x(:,1)-r_eq)/1000;
E = 0.5*x(:,4).^2 - planet.mu./x(:,1); %Actual energy
En = (E-E(1))/(E(end)-E(1)); %Normalized Energy, 0 at entry and 1 at deployment

[theta_t,phi_t] = FinalLatLon(x(1,2),x(1,3),x(1,6),DR_t,CR_t);
output.target.lat = phi_t;
output.target.lon = theta_t;
output.target.DR = DR_t;
output.target.CR = CR_t;
output.energy = E;
output.energy_norm = En;

% Altitude vs Velocity
figure
ParachuteDeploymentConstraints(true);
hold on
plot(x(:,4),hkm)
xlabel('Velocity (m/s)')
ylabel('Altitude (km)')
title(['Final altitude: ',num2str(hkm(end)), ' km'])

%DR v CR
[DR,CR] = Range(x(1,2),x(1,3),x(1,6),x(:,2),x(:,3));
output.DR = DR;
output.CR = CR;
figure
plot(CR,DR)
ylabel('Downrange (km)')
xlabel('Crossrange (km)')
title(['Downrange: ',num2str(DR(end)),' km, Crossrange: ',num2str(CR(end)), ' km'])

%Lat/Lon
figure
plot(x(1,3)/dtr,x(1,2)/dtr,'ko')
hold on
plot(phi_t/dtr,theta_t/dtr,'kx')
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
plot(En,x(:,5)/dtr)
xlabel('Normalized Energy')
ylabel('FPA (deg)')

%Heading angle
figure
subplot 211
plot(t,x(:,6)/dtr)
xlabel('Time (s)')
ylabel('Heading (deg)')
subplot 212
plot(En,x(:,6)/dtr)
xlabel('Normalized Energy')
ylabel('Heading (deg)')
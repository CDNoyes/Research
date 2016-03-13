%ENTRYPLOTS Plot some standard graphs useful in examining entry
%trajectories.

function EntryPlots(trajectorySummary)
dtr = pi/180;
planet = Mars();
ts = trajectorySummary;
t = ts.time;
x = ts.state;
sig = ts.sigma;
u = ts.control;
r_eq = planet.radiusEquatorial;
hkm = (x(:,1)-r_eq)/1000;

lineWidth = 2;
markerSize = 10;
fontSize = 12;
fontWeight = 'bold';
fontColor = 'k';

% Altitude vs Velocity
figure
ParachuteDeploymentConstraints(true);
hold on
plot(x(:,4),hkm, 'LineWidth',lineWidth)
xlabel('Velocity (m/s)','FontSize',fontSize,'FontWeight',fontWeight,'Color',fontColor)
ylabel('Altitude (km)','FontSize',fontSize,'FontWeight',fontWeight,'Color',fontColor)
title(['Final altitude: ',num2str(hkm(end)), ' km'],'FontSize',fontSize,'FontWeight',fontWeight,'Color',fontColor)

%DR v CR
figure
plot(ts.CR,ts.DR, 'LineWidth',lineWidth)
ylabel('Downrange (km)','FontSize',fontSize,'FontWeight',fontWeight,'Color',fontColor)
xlabel('Crossrange (km)','FontSize',fontSize,'FontWeight',fontWeight,'Color',fontColor)
title(['Downrange: ',num2str(ts.DR(end)),' km, Crossrange: ',num2str(ts.CR(end)), ' km'],'FontSize',fontSize,'FontWeight',fontWeight,'Color',fontColor)

%Lat/Lon
figure
plot(x(1,3)/dtr,x(1,2)/dtr,'ko','MarkerSize',markerSize)
hold on
plot(ts.target.lat/dtr,ts.target.lon/dtr,'kx','MarkerSize',markerSize)
plot(x(:,3)/dtr,x(:,2)/dtr, 'LineWidth',lineWidth)
xlabel('Latitude (deg)','FontSize',fontSize,'FontWeight',fontWeight,'Color',fontColor)
ylabel('Longitude (deg)','FontSize',fontSize,'FontWeight',fontWeight,'Color',fontColor)
title('Ground Track', 'FontSize',fontSize,'FontWeight',fontWeight,'Color',fontColor)
legend('Initial point','Target','Location','Best')

%Flight path angle
figure
subplot 211
plot(t,x(:,5)/dtr, 'LineWidth',lineWidth)
axis([0,max(t),min(x(:,5)/dtr)*(1-sign(min(x(:,5)))*.1),max(x(:,5)/dtr)*(1+sign(max(x(:,5)))*.1)])
xlabel('Time (s)','FontSize',fontSize,'FontWeight',fontWeight,'Color',fontColor)
ylabel('FPA (deg)','FontSize',fontSize,'FontWeight',fontWeight,'Color',fontColor)
subplot 212
plot(ts.energy_norm,x(:,5)/dtr, 'LineWidth',lineWidth)
axis([0,1,min(x(:,5)/dtr)*(1-sign(min(x(:,5)))*.1),max(x(:,5)/dtr)*(1+sign(max(x(:,5)))*.1)])
xlabel('Normalized Energy','FontSize',fontSize,'FontWeight',fontWeight,'Color',fontColor)
ylabel('FPA (deg)','FontSize',fontSize,'FontWeight',fontWeight,'Color',fontColor)

%Heading angle
figure
subplot 211
plot(t,x(:,6)/dtr, 'LineWidth',lineWidth)
axis([0,max(t),min(x(:,6)/dtr)*(1-sign(min(x(:,6)))*.1),max(x(:,6)/dtr)*(1+sign(max(x(:,6)))*.1)])
xlabel('Time (s)','FontSize',fontSize,'FontWeight',fontWeight,'Color',fontColor)
ylabel('Heading (deg)','FontSize',fontSize,'FontWeight',fontWeight,'Color',fontColor)
subplot 212
plot(ts.energy_norm,x(:,6)/dtr, 'LineWidth',lineWidth)
axis([0,1,min(x(:,6)/dtr)*(1-sign(min(x(:,6)))*.1),max(x(:,6)/dtr)*(1+sign(max(x(:,6)))*.1)])
xlabel('Normalized Energy','FontSize',fontSize,'FontWeight',fontWeight,'Color',fontColor)
ylabel('Heading (deg)','FontSize',fontSize,'FontWeight',fontWeight,'Color',fontColor)

%CONTROL
figure
subplot 211
plot(t,sig/dtr, 'LineWidth',lineWidth)
axis([0,max(t),-90,90])
xlabel('Time (s)','FontSize',fontSize,'FontWeight',fontWeight,'Color',fontColor)
ylabel('Bank Angle (deg)','FontSize',fontSize,'FontWeight',fontWeight,'Color',fontColor)
subplot 212
plot(ts.energy_norm,sig/dtr, 'LineWidth',lineWidth)
axis([0,1,-90,90])
xlabel('Normalized Energy','FontSize',fontSize,'FontWeight',fontWeight,'Color',fontColor)
ylabel('Bank Angle (deg)','FontSize',fontSize,'FontWeight',fontWeight,'Color',fontColor)

% Drag dynamics, for analysis
if true
    figure
    subplot 211
    title('Drag Dynamics')
    plot(ts.time,[ts.a; ts.b ],'--','LineWidth',2)
    hold all
    plot(ts.time,ts.a + ts.b .* ts.control,'LineWidth',2)
    legend('a','b','a+bu')
    axis([0,max(t),-2,1])
    xlabel('Time (s)','FontSize',fontSize,'FontWeight',fontWeight,'Color',fontColor)
    
    subplot 212
    plot(ts.energy_norm,[ts.a; ts.b],'LineWidth',2)
    hold all
    plot(ts.energy_norm,ts.a + ts.b .* ts.control,'LineWidth',2)
    axis([0,1,-2,1])
    xlabel('Normalized Energy','FontSize',fontSize,'FontWeight',fontWeight,'Color',fontColor)
    
end
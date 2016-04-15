%ENTRYPLOTS Plot some standard graphs useful in examining entry
%trajectories.

function EntryPlots(trajectorySummary)
% close all
dtr = pi/180;
planet = Mars();
ts = trajectorySummary;
t = ts.time;
x = ts.state;
sig = ts.sigma;
u = ts.control;
r_eq = planet.radiusEquatorial;
hkm = (x(:,1)-r_eq)/1000;
e = ts.energy_norm;

lineWidth = 2;
markerSize = 10;
fontSize = 12;
fontWeight = 'bold';
fontColor = 'k';

% Altitude vs Velocity
n = 1;
figure(n)
grid on
box on
set(gcf,'name','Altitude vs Velocity', 'numbertitle','off','WindowStyle','docked')
ParachuteDeploymentConstraints(true);
hold on
plot(x(:,4),hkm, 'LineWidth',lineWidth)
xlabel('Velocity (m/s)','FontSize',fontSize,'FontWeight',fontWeight,'Color',fontColor)
ylabel('Altitude (km)','FontSize',fontSize,'FontWeight',fontWeight,'Color',fontColor)
title(['Final altitude: ',num2str(hkm(end)), ' km'],'FontSize',fontSize,'FontWeight',fontWeight,'Color',fontColor)

%DR v CR
n = n+1;
figure(n)
plot(ts.CR,ts.DR, 'LineWidth',lineWidth)
ylabel('Downrange (km)','FontSize',fontSize,'FontWeight',fontWeight,'Color',fontColor)
xlabel('Crossrange (km)','FontSize',fontSize,'FontWeight',fontWeight,'Color',fontColor)
title(['Downrange: ',num2str(ts.DR(end)),' km, Crossrange: ',num2str(ts.CR(end)), ' km'],'FontSize',fontSize,'FontWeight',fontWeight,'Color',fontColor)
grid on
box on
set(gcf,'name','Downrange vs Crossrange', 'numbertitle','off','WindowStyle','docked')

%Lat/Lon
n = n+1;
figure(n)
plot(x(1,3)/dtr,x(1,2)/dtr,'ko','MarkerSize',markerSize)
hold on
plot(ts.target.lat/dtr,ts.target.lon/dtr,'kx','MarkerSize',markerSize)
plot(x(:,3)/dtr,x(:,2)/dtr, 'LineWidth',lineWidth)
xlabel('Latitude (deg)','FontSize',fontSize,'FontWeight',fontWeight,'Color',fontColor)
ylabel('Longitude (deg)','FontSize',fontSize,'FontWeight',fontWeight,'Color',fontColor)
title('Ground Track', 'FontSize',fontSize,'FontWeight',fontWeight,'Color',fontColor)
legend('Initial point','Target','Location','Best')
grid on
box on
set(gcf,'name','Groundtrack', 'numbertitle','off','WindowStyle','docked')

%Flight path angle
n = n+1;
figure(n)
subplot 211
plot(t,x(:,5)/dtr, 'LineWidth',lineWidth)
grid on
axis([0,max(t),min(x(:,5)/dtr)*(1-sign(min(x(:,5)))*.1),max(x(:,5)/dtr)*(1+sign(max(x(:,5)))*.1)])
xlabel('Time (s)','FontSize',fontSize,'FontWeight',fontWeight,'Color',fontColor)
ylabel('FPA (deg)','FontSize',fontSize,'FontWeight',fontWeight,'Color',fontColor)
subplot 212
plot(ts.energy_norm,x(:,5)/dtr, 'LineWidth',lineWidth)
axis([0,1,min(x(:,5)/dtr)*(1-sign(min(x(:,5)))*.1),max(x(:,5)/dtr)*(1+sign(max(x(:,5)))*.1)])
xlabel('Normalized Energy','FontSize',fontSize,'FontWeight',fontWeight,'Color',fontColor)
ylabel('FPA (deg)','FontSize',fontSize,'FontWeight',fontWeight,'Color',fontColor)
grid on
box on
set(gcf,'name','FPA', 'numbertitle','off','WindowStyle','docked')

%Heading angle
n = n+1;
figure(n)
subplot 211
plot(t,x(:,6)/dtr, 'LineWidth',lineWidth)
grid on
axis([0,max(t),min(x(:,6)/dtr)*(1-sign(min(x(:,6)))*.1),max(x(:,6)/dtr)*(1+sign(max(x(:,6)))*.1)])
xlabel('Time (s)','FontSize',fontSize,'FontWeight',fontWeight,'Color',fontColor)
ylabel('Heading (deg)','FontSize',fontSize,'FontWeight',fontWeight,'Color',fontColor)
subplot 212
plot(ts.energy_norm,x(:,6)/dtr, 'LineWidth',lineWidth)
axis([0,1,min(x(:,6)/dtr)*(1-sign(min(x(:,6)))*.1),max(x(:,6)/dtr)*(1+sign(max(x(:,6)))*.1)])
xlabel('Normalized Energy','FontSize',fontSize,'FontWeight',fontWeight,'Color',fontColor)
ylabel('Heading (deg)','FontSize',fontSize,'FontWeight',fontWeight,'Color',fontColor)
grid on
box on
set(gcf,'name','Heading', 'numbertitle','off','WindowStyle','docked')

%CONTROL
n = n+1;
figure(n)
subplot 211
plot(t,Saturate(sig/dtr,-90,90), 'LineWidth',lineWidth)
grid on
axis([0,max(t),-95,95])
xlabel('Time (s)','FontSize',fontSize,'FontWeight',fontWeight,'Color',fontColor)
ylabel('Bank Angle (deg)','FontSize',fontSize,'FontWeight',fontWeight,'Color',fontColor)
subplot 212
plot(ts.energy_norm,Saturate(sig/dtr,-90,90), 'LineWidth',lineWidth)
axis([0,1,-95,95])
xlabel('Normalized Energy','FontSize',fontSize,'FontWeight',fontWeight,'Color',fontColor)
ylabel('Bank Angle (deg)','FontSize',fontSize,'FontWeight',fontWeight,'Color',fontColor)
grid on
box on
set(gcf,'name','Bank Profile', 'numbertitle','off','WindowStyle','docked')



if isfield(ts,'observer') && ts.observer
    n = n+1;
    figure(n)
    plot(ts.energy_norm,ts.state(:,7),'LineWidth',lineWidth)
    hold all
    plot(ts.energy_norm,ts.D,'LineWidth',lineWidth)
    plot(ts.energy_norm,ts.D'./cos(x(:,5)), 'LineWidth',lineWidth)
    
    legend('Observed Drag','Modeled Drag', 'Modeled Drag/cos(\gamma)')
    xlabel('Normalized Energy (-)','FontSize',fontSize,'FontWeight',fontWeight,'Color',fontColor)
    ylabel('Drag (m/s^2)','FontSize',fontSize,'FontWeight',fontWeight,'Color',fontColor)
    set(gcf,'name','Drag Profile', 'numbertitle','off','WindowStyle','docked')
    grid on
    
    n = n+1;
    figure(n)
    plot(ts.energy_norm,ts.state(:,9),'LineWidth',lineWidth)
    xlabel('Normalized Energy (-)','FontSize',fontSize,'FontWeight',fontWeight,'Color',fontColor)
    ylabel('Disturbance Estimate (m/s^4)','FontSize',fontSize,'FontWeight',fontWeight,'Color',fontColor)
    grid on
    set(gcf,'name','Disturbance Estimate', 'numbertitle','off','WindowStyle','docked')
end

%Experimental - Reference
% n = n+1;
% figure(n)
% subplot 211
% plot(ts.time,ts.D'./cos(x(:,5)), 'LineWidth',lineWidth)
% grid on
% title('Drag divided by cos(fpa)')
% axis([0,max(t),-10,90])
% xlabel('Time (s)','FontSize',fontSize,'FontWeight',fontWeight,'Color',fontColor)
% ylabel(' (m/s^2)','FontSize',fontSize,'FontWeight',fontWeight,'Color',fontColor)
% subplot 212
% plot(ts.energy_norm,ts.D'./cos(x(:,5)), 'LineWidth',lineWidth)
% axis([0,1,-10,90])
% xlabel('Normalized Energy','FontSize',fontSize,'FontWeight',fontWeight,'Color',fontColor)
% ylabel('m/s^2','FontSize',fontSize,'FontWeight',fontWeight,'Color',fontColor)
% grid on
% box on
% set(gcf,'name','Reference Profile', 'numbertitle','off','WindowStyle','docked')


% Various experimental crossrange plans
if true
    n = n+1;
    figure(n)
    plot(ts.L./ts.D,abs(ts.CR),'LineWidth',lineWidth)
    set(gcf,'name','Crossrange vs L/D', 'numbertitle','off','WindowStyle','docked')
    
    n = n+1;
    figure(n)
    plot(ts.D, abs(ts.CR),'LineWidth',lineWidth)
    set(gcf,'name','Crossrange vs Drag', 'numbertitle','off','WindowStyle','docked')
    grid on
    
    n = n+1;
    figure(n)
    plot(ts.L'.*cos(x(:,5)), abs(ts.CR),'LineWidth',lineWidth)
    set(gcf,'name','Crossrange vs Lcos\gamma', 'numbertitle','off','WindowStyle','docked')
    grid on
    
    n = n+1;
    figure(n)
    plot(ts.state(:,4), abs(ts.CR),'LineWidth',lineWidth)
    set(gcf,'name','Crossrange vs Velocity', 'numbertitle','off','WindowStyle','docked')
    grid on
end

% Drag dynamics, for analysis
if false
    n = n+1;
    figure(n)
    plot(ts.energy_norm,ts.D, 'LineWidth',lineWidth)
    axis([0,1,0,1.1*max(ts.D)])
    xlabel('Normalized Energy','FontSize',fontSize,'FontWeight',fontWeight,'Color',fontColor)
    ylabel('Drag (m/s^2)','FontSize',fontSize,'FontWeight',fontWeight,'Color',fontColor)
    
    n = n+1;
    figure(n)
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
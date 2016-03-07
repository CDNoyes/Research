%Runs the van der Pol example 
clean
x0 = [.2;.1];
z0 = 0;

[t,x] = ode45(@TerminalSMCDynamics,[0,50],[x0;z0]);
for i = 1:length(t)
    [~,u(i),d(i)] = TerminalSMCDynamics(t(i),x(i,:));
end
figure
plot(t,x(:,1),t,2*sin(0.2*t))
% hold all
% plot(t,x(:,2),t,.4*cos(0.2*t))

figure
plot(t,x(:,1)-2*sin(0.2*t))

figure
plot(t,u)

figure
plot(t,d)
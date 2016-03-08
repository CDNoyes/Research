function TimeVaryingSMC()

x0 = [0.5;0.4];
e0 = [x0(1),x0(2)-0.4];
[t,x] = ode45(@dynamics,[0,100],x0,[],e0);
for i = 1:length(t)
    [~,u(i),s(i),d(i),T] = dynamics(t(i),x(i,:)',e0);
end
iT = find(t>T,1);

figure
plot(x(:,1)-2*sin(0.2*t),x(:,2)-0.2*2*cos(0.2*t))
hold on
plot(x(iT,1)-2*sin(0.2*t(iT)),x(iT,2)-0.2*2*cos(0.2*t(iT)),'ro')
plot(0,0,'kx',e0(1),e0(2),'ko')
legend('Trajectory',['Time invariance (t=',num2str(t(iT)),' s)'],'Target','IC')
title('Error trajectory')

figure
plot(t,[2*sin(0.2*t),0.2*2*cos(0.2*t)])
hold all
plot(t,d)
plot(t,x)
legend('Ref','Ref rate','Disturbance')

figure
plot(t,u)
title(['u_0 = ',num2str(u(1)),', u_{max} = ',num2str(max(abs(u)))])
ylabel('u')

figure
plot(t,s)
title('Sliding Surface')
axis([0,t(end),-1,1])
end

function [dx,u,s,d,T] = dynamics(t,x,e0)


%System format
f = -2*x(1)+3*(1-x(1)^2)*x(2);
b = 1;

% Reference Trajectory
yd = 2*sin(0.2*t);
yd_dot = 0.2*2*cos(0.2*t);
yd_ddot = -0.2^2*2*sin(0.2*t);

% Disturbance
d = 2*sin(0.1*pi*t)+3*sin(0.2*sqrt(t+1));
% d = 0;
dmax = d; %Estimated bound of the disturbance

% Error states
e = x-[yd;yd_dot];

% Sliding Surface Design Parameters
[T,A,B,~,beta] = DesignSM(e0(1),8.8);
if t == 0
    disp(beta)
end
c = [beta,1];

%Sliding Surface Computation
s = c*e + (A*t+B)*(t<=T); 

% Control
u0 = yd_ddot-f - dmax*Saturate(s/.01,-1,1);
U = 10;
u = Saturate(1/b*u0 + 1/b*(-beta*e(2)-A*(t<=T)),-U,U);

dx = [  x(2)
        f + b*u + d];
    
end
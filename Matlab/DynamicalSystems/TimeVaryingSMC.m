% Reference: "Time Varying Sliding Modes for Second Order Systems"


function TimeVaryingSMC()

%Initial States
x0 = [0.9;0.4];
x0_hat = [0.1;0.4];
e0 = [x0(1),x0(2)-0.4];
d0 = 0.6;

%Design controller, and observer
% Sliding Surface Design Parameters
% [S,controller,T] = DesignSM(e0(1),1,'Terminal');
[S,controller,T] = DesignSM(e0(1),1,'Velocity');
% [S,controller,T] = DesignSM(e0(1),1,'Acceleration');
[alpha, k] = computeSMOGains(400);

%Integrate the controlled trajectory
[t,x] = ode45(@dynamics,[0,100],[x0;x0_hat;d0],[],S,controller,k,alpha);
for i = 1:length(t)
    [~,u(i),s(i),d(i),d_hat(i)] = dynamics(t(i),x(i,:)',S,controller,k,alpha);
end
iT = find(t>T,1);

%Plot the results
figure
plot(x(:,1)-2*sin(0.2*t),x(:,2)-0.2*2*cos(0.2*t))
hold on
plot(x(iT,1)-2*sin(0.2*t(iT)),x(iT,2)-0.2*2*cos(0.2*t(iT)),'ro')
plot(0,0,'kx',e0(1),e0(2),'ko')
legend('Trajectory',['Time Invariance (t = ',num2str(T),' s)'],'Target','IC')
title('Error trajectory')

figure
plot(t,[2*sin(0.2*t),0.2*2*cos(0.2*t)])
hold all
plot(t,d,'--')
plot(t,x(:,1:2))
legend('Ref','Ref rate','Disturbance')

figure
plot(t,u)
title(['u_0 = ',num2str(u(1)),', u_{max} = ',num2str(max(abs(u)))])
ylabel('u')

figure
plot(t,s)
title('Sliding Surface')
axis([0,t(end),-1,1])

% figure
% plot(t,d)
% hold all
% plot(t,d_hat,t,x(:,6))
% title('Disturbance Estimate History')

figure
semilogy(t,abs(x(:,1:2)-x(:,3:4)))
hold all
semilogy(t,abs(d-d_hat))
title('Estimation Error History')
legend('x_1','x_2','d')
end

function [dX,u,s,d,d_hat] = dynamics(t,X,S,controller,k,alpha)
%Parse the extended state vector
x = X(1:2);
x_hat = X(3:4);
d_hat = X(5);

%System format
f = -2*x(1)+3*(1-x(1)^2)*x(2);
b = 1;

% Reference Trajectory
yd = 2*sin(0.2*t);
yd_dot = 0.2*2*cos(0.2*t);
yd_ddot = -0.2^2*2*sin(0.2*t);

% Disturbance
d = 2*sin(0.1*pi*t)+3*sin(0.2*sqrt(t+1));
% d = d/10;
% d = 0;
dmax = 5*(t<.1) + abs(d_hat)*(t>=.1); % Estimated bound of the disturbance

% Error states
e = x-[yd;yd_dot];

%Sliding Surface Computation
s = S(e,t);

% Control
u0 = yd_ddot-f - dmax*Saturate(s/.01,-1,1);
u = controller(b,u0,e,t);
U = 100;
u = Saturate(u,-U,U);


dx = [  x(2)
        f + b*u + d]; 

    % This is the disturbance observer term
dx_hat = SMO([x_hat;d_hat],x(1),f,b,u,k,alpha);

dX = [dx;dx_hat];

end

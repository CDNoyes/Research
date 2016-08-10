%Convert constraints on dynamic pressure and Mach number into altitude and
%velocity constraints
function [v1,v2,v3,h_min,m] = ParachuteDeploymentConstraints(PLOT)
if nargin == 0
    PLOT = false;
end
q1 = 300;
q2 = 800;

M1 = 1.4;
M2 = 2.2;

h_min = 6e3; %min altitude, meters, artificial limit based on mission req.

h = linspace(h_min, 18e3,500);

for i = 1:length(h)
    [rho,a] = MarsAtmosphericDensity(h(i),0); 
    Vmax(i) = min(a*M2,sqrt(2*q2/rho));                               
    Vmin(i) = max(a*M1,sqrt(2*q1/rho));   
    if Vmin(i) >= Vmax(i)
        j = i;
        break
    end
end
[Max,Index] = max(Vmax);

x = linspace(Vmin(1),Vmax(1));
y = ones(1,length(x))*h_min/1000;
h = h(1:j);
if PLOT
    [lineSpecs,textSpecs,figSpecs] = PlotSpecs();
plot(Vmax,h/1000,Vmin,h/1000,x,y,lineSpecs{:})
xlabel('Velocity (m/s)',textSpecs{:})
ylabel('Altitude (km)',textSpecs{:})

end
%find the equation of the line
m = (h(Index)-h_min)/(Max-Vmax(1)); %slope
% vtest = linspace(Vmin(1),Max);
% for i = 1:length(vtest)
%     if vtest(i) <= Vmax(1)
%         htest(i) = h_min;
%     else
%         htest(i) = h_min+m*(vtest(i)-Vmax(1));
%     end
% end
% %test plot
% plot(vtest,htest/1000)
v1 = Vmin(1);
v2 = Vmax(1);
v3 = Max;
end
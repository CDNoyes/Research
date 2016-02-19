%TRAJECTORY Creates a structure with all the data possibly
%required for trajectory tracking computations.

function output = TrajectorySummary(t,x,DR_t,CR_t)

tf = EntryAnalysis(t,x,DR_t,CR_t);
t = t(t<=tf);
x = x(1:length(t),:);

% dtr = pi/180;
planet = Mars();
vm = VehicleModel();
% r_eq = planet.radiusEquatorial;
% hkm = (x(:,1)-r_eq)/1000;
E = 0.5*x(:,4).^2 - planet.mu./x(:,1); %Actual energy
En = (E-E(1))/(E(end)-E(1)); %Normalized Energy, 0 at entry and 1 at deployment

[theta_t,phi_t] = FinalLatLon(x(1,2),x(1,3),x(1,6),DR_t,CR_t);
output.target.lat = phi_t;
output.target.lon = theta_t;
output.target.DR = DR_t;
output.target.CR = CR_t;
output.energy = E;
output.energy_norm = En;
output.time = t;
output.state = x;
[DR,CR] = Range(x(1,2),x(1,3),x(1,6),x(:,2),x(:,3));
output.DR = DR;
output.CR = CR;

for i = 1:length(t)
    [output.g(i),output.L(i),output.D(i),~,output.M(i),output.a(i),output.rho(i)] = EntryForces(x(i,:),planet,vm);
end

end
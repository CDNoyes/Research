function EntryAnalysis(t,x)

dtr = pi/180;
planet = Mars();
r_eq = planet.radiusEquatorial;
hkm = (x(:,1)-r_eq)/1000;
E = 0.5*x(:,4).^2 - planet.mu./x(:,1); %Actual energy
En = (E-E(1))/(E(end)-E(1)); %Normalized Energy, 0 at entry and 1 at deployment
[DR,CR] = Range(x(1,2),x(1,3),x(1,6),x(:,2),x(:,3));
tf = findTrajLength(t,x);


disp(['Final altitude: ',num2str(hkm(end)), ' km'])
disp(['Final FPA: ',num2str(x(end,5)/dtr), ' deg'])
disp(['Final Downrange: ',num2str(DR(end)),' km'])
disp(['Final Crossrange: ',num2str(CR(end)), ' km'])
disp(['Trajectory Duration: ',num2str(tf), ' s'])

end

function tf = findTrajLength(t,x)

d = diff(x(:,1));
idx = find(d==0,1);
tf = t(idx);

end
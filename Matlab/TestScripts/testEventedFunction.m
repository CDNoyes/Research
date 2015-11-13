clear
f = @(t,x) [-x(1)-1; -x(2)/(x(1)+2)];
x0 = [1;1];

[t,x] = ode45(f,[0,5],x0);

%with "event"
opt = odeset('AbsTol',1e-12,'RelTol',1e-12);
[te,xe] = ode45( @(t,x) eventedFunction(t,x,f,'t > 2 || x(1) < -0.1'), [0,5], x0,opt);

id = find(diff(xe)==0,1);

plot(t,x)
hold all
plot(te(1:id),xe(1:id,:),'r--')
plot(te(id),xe(id,:),'ko')
line(ones(1,100)*te(id),linspace(min(xe(id,:)),max(xe(id,:))),'Color','k','LineStyle','--')
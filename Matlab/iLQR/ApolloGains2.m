function K = ApolloGains2(t, v, x, u, L, D)

h = x(1,:);
fpa = x(2,:);
s = x(3,:);

r = h + 3396.2e3;
cg = cos(fpa);
sg = sin(fpa);
akm = 1000;
g = 3.71;
hs = 9354.5;
liftv = L.*u;

c1 = liftv./v.^2 + cg./r + g*cg./v.^2;
c2 = (v./r-g./v).*sg;

N = length(v);
L = zeros(5,N); 
L(:,end) = [1; 0; 0; -1/tan(fpa(end)); 0]; 

dt = diff(t);

for i = N:-1:2
    l = L(:,i);
    dl1 = 0;
    dl2 = -cg(i)*l(1) + 2*D(i)/v(i)*l(2) - c1(i)*l(3)-sg(i)*l(4);
    dl3 = v(i)*sg(i)*l(1) + g*cg(i)*l(2) + c2(i)*l(3) - v(i)*cg(i)*l(4);
    dl4 = (-D(i)/hs-2*g/r(i))*l(2) + l(3)*(liftv(i)/v(i)/hs + v(i)*cg(i)/r(i).^2 - cg(i)*2*g/r(i)/v(i));
    dl5 = -D(i)/v(i)*l(3);
    
    dL = [dl1;dl2;dl3;dl4;dl5];
    L(:,i-1) =  l - dL*dt(i-1);   
end

% Reorder for [d,s,f]
K = -[L(4,:).*-hs./D; L(1,:); L(3,:)]./repmat(L(5,:),3,1); % need to convert altitude gain to drag
% s gain is negated because we're using range flown instead of range to go 

if 0
    for i = 1:5
        figure
        plot(v, L(i,:))
    end
end
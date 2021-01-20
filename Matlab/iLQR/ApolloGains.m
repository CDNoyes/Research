function K = ApolloGains(t, v, x, u, L, D, A, B)

h = x(1,:);
fpa = x(2,:);
s = x(3,:);

% A and B are in the order [v, h, fpa, s]
N = length(v);
F = zeros(4,N); % v h fpa s u
F(:,end) = [0; -1/tan(fpa(end)); 0; 1]; 
Fu(N) = 0;

dt = diff(t);

for i = N:-1:2
    Acl = A(:,:,i);
    F(:,i-1) =  F(:,i) + Acl.'*F(:,i)*dt(i-1);
    Fu(i-1) =  Fu(i) + L(i)/v(i)*F(3,i)*dt(i-1); % u as the control
    
end
Fu(Fu==0) = 0.1; %min(abs(Fu(Fu~=0)));

% Reorder for [d,s,f]

K = -[F(2,:); -F(4,:); F(3,:)]./repmat(Fu,3,1); % need to convert altitude gain to drag outside
% K = -[F(2,:); F(4,:); F(3,:); Fu]; % need to convert altitude gain to drag outside
% s gain will also need to be negated because we're using range flown instead of range to go 

if 0
    for i = 1:3
        figure
        plot(v, F(i,:))
    end
    figure
    plot(v, Fu)
end
%% Variable mu

close all
clc
clear

[lineSpecs,textSpecs,figSpecs] = PlotSpecs();


% generate actual integration data to compare to
[F,J] = VDP();
f = @(t,x,u,mu) [F(t,x,u,mu); -x(3)*trace(J(x,u,mu)*[eye(2);zeros(2)]) ];
opt = odeset('AbsTol',1e-8,'RelTol',1e-8);

muNom = 0.5;
uniform = 1;
if uniform
    mus = linspace(0.25,0.75,40);
    p0 = ones(size(mus))/(mus(end)-mus(1));
    text = ['U','[',num2str(mus(1)),', ',num2str(mus(end)),']'];
else %normal
%     mus = sort(muNom+.09*randn(1,100)); %Sample from distribution
    mus = linspace(muNom-3.1*0.09,muNom+3.1*0.09,40);
    p0 = pdf('Normal',mus,muNom,0.09);
    text = ['N(',num2str(muNom),', ',num2str(0.09),')'];
    figure(10)
    plot(mus,p0,'*')
end

x0 = [0;-.1];
tf = 8;
[t,x] = ode45(f,[0,tf],[x0;0],opt,@(t) 0, muNom);
h0 = plot(x(:,1),x(:,2),'k',lineSpecs{:});
n = length(t);
% h = logspace(-1,-3,n);
% t = (max(h)-h)*tf/(h(1)-h(end));

X = zeros(n,2);
E = zeros(n,2);
P = zeros(n,1);
[Xmat{1:2},State{1:3}] = deal(zeros(n,length(mus)));
figure(1)
title(['\mu \in', text],textSpecs{:})
hold all

for i = 1:length(mus)
    [T,x] = ode45(f,t,[x0;p0(i)],opt,@(t) 0, mus(i));
    if i == 1 || i==length(mus)
        plot(x(:,1),x(:,2),'g*')
    else
        plot(x(:,1),x(:,2),'c')
    end
    plot(x(end,1),x(end,2),'m*')

    X = X + x(:,1:2); %Average trajectory
    E = E + [x(:,1).*x(:,3),x(:,2).*x(:,3)];
    P = P + x(:,3);

    State{1}(:,i) = x(:,1);
    State{2}(:,i) = x(:,2);
    State{3}(:,i) = x(:,3);

    xf(i,:) = x(end,:);
end


h1 = plot(X(:,1)/length(mus),X(:,2)/length(mus),'b',lineSpecs{:});
h2 = plot(E(:,1)./P,E(:,2)./P,'r',lineSpecs{:});
% h3 = plot(E2(:,1),E2(:,2),'g',lineSpecs{:});
xlabel('x_1',textSpecs{:})
ylabel('x_2',textSpecs{:})
legend([h0 h1, h2],{'f(E[\mu])','E[f(\mu)] via Sampling','E[f(\mu)] via PF'})

figure(3)
surf(State{1},State{2},repmat(mus,n,1),State{3})
xlabel('x_1',textSpecs{:})
ylabel('x_2',textSpecs{:})
zlabel('\mu',textSpecs{:})

figure(4)
surf(State{1},State{2},State{3},repmat(mus,n,1))
xlabel('x_1',textSpecs{:})
ylabel('x_2',textSpecs{:})
zlabel('p',textSpecs{:})


for i = [1:10:length(T),length(T)]
    figure(5)
    hold all
    plot3(State{1}(i,:), T(i)*ones(1,length(mus)), State{3}(i,:),'*')
    figure(6)
    hold all
    plot3(State{2}(i,:), T(i)*ones(1,length(mus)), State{3}(i,:),'*')
end


return

%%

x0_nom = [1;-1];
n = 500;
%sigma = [.5 .1;.1 .5];
sigma = diag([1 1]);
%sigma = [.9 .1; .1 0.1];
x0_rand = mvnrnd(x0_nom',sigma,n);
p0 = mvnpdf(x0_rand, x0_nom',sigma);

scatter3(x0_rand(:,1),x0_rand(:,2),p0,[],p0)


x0_rand = [x0_rand,p0];

for i = 1:1:n
    [T,x] = ode45(f,t,x0_rand(i,:),[],@(t) 0, 0.5);
    State{1}(:,i) = x(:,1);
    State{2}(:,i) = x(:,2);
    State{3}(:,i) = x(:,3);
end

figure(5)
surf(State{1},State{2},State{3})
xlabel('x_1',textSpecs{:})
ylabel('x_2',textSpecs{:})
zlabel('p',textSpecs{:})


%%
e = [0.0188 0.0123
    .0109 .0284
    0.0269 0.029
    .0103 0.0449
    .0309 0.0172
    0.0255 0.0171]*6378*pi/180; 
for i = 1:6
    norm(e(i,:))
end
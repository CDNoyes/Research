%% Uncertain mu

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
    mus = linspace(0.25,0.75,300);
    musMC = mus;
    p0 = ones(size(mus))/(mus(end)-mus(1));
    text = ['U','[',num2str(mus(1)),', ',num2str(mus(end)),']'];
    boundLabel = 'True uniform bounds';
else %normal
    stdDev = 0.08;
    musMC = sort(muNom+stdDev*randn(1,50)); %Sample from distribution for MC
    mus = linspace(muNom-3*stdDev,muNom+3*stdDev,10); %PF has to be nicely spaced over the whole domain
    p0 = pdf('Normal',mus,muNom,stdDev);
    text = ['N(',num2str(muNom),', ',num2str(stdDev),')'];
    boundLabel = 'True 3-sigma bounds';
    figure(10)
    plot(musMC,pdf('Normal',musMC,muNom,stdDev),'*')
end

x0 = [0;-.1];
tf = 8;
[t,x] = ode45(f,[0,tf],[x0;0],opt,@(t) 0, muNom);
figure(1)
h0 = plot(x(:,1),x(:,2),'k',lineSpecs{:});
n = length(t);


X = zeros(n,2);
[State{1:3}] = deal(zeros(n,length(mus)));
figure(1)
title(['\mu \in', text],textSpecs{:})
hold all

for i = 1:length(mus)
    [T,x] = ode45(f,t,[x0;p0(i)],opt,@(t) 0, mus(i));
    if i == 1 || i==length(mus)
        hb = plot(x(:,1),x(:,2),'g*');
    else
        %         plot(x(:,1),x(:,2),'c')
    end
    plot(x(end,1),x(end,2),'m*')
    if uniform
        X = X + x(:,1:2); %Average trajectory, not a proper approximation for Gaussian mu
        xf(i,:) = x(end,:);
        
    end
    State{1}(:,i) = x(:,1);
    State{2}(:,i) = x(:,2);
    State{3}(:,i) = x(:,3);
    
end
Ptotal = trapz(mus,State{3},2);
State{3} = State{3}./repmat(Ptotal,1,length(mus));
% E = [trapz(mus,State{1}.*State{3},2), trapz(mus,State{2}.*State{3},2)];
% V = [trapz(mus,(State{1}-repmat(E(:,1),1,length(mus))).^2.*State{3},2), trapz(mus,(State{2}-repmat(E(:,2),1,length(mus))).^2.*State{3},2)];
% V(:,3) = trapz(mus,(State{1}-repmat(E(:,1),1,length(mus))).*(State{2}-repmat(E(:,2),1,length(mus))).*State{3},2);
%
E = [trapz(mus,State{1}.*repmat(p0,length(t),1),2), trapz(mus,State{2}.*repmat(p0,length(t),1),2)];
V = [trapz(mus,(State{1}-repmat(E(:,1),1,length(mus))).^2.*repmat(p0,length(t),1),2), trapz(mus,(State{2}-repmat(E(:,2),1,length(mus))).^2.*repmat(p0,length(t),1),2)];
V(:,3) = trapz(mus,(State{1}-repmat(E(:,1),1,length(mus))).*(State{2}-repmat(E(:,2),1,length(mus))).*repmat(p0,length(t),1),2);

if ~uniform
    for i = 1:length(musMC)
        [T,x] = ode45(F,t,x0,opt,@(t) 0, musMC(i));
        X = X + x(:,1:2); %Average trajectory
        xf(i,:) = x(end,:);
        
    end
    X = X/length(musMC);
    for i = 1:10:length(T)
        P = diag(V(i,1:2)) + V(i,3)*(ones(2)-eye(2));
        DrawNormalEllipse(E(i,:)+.0*T(i)*[State{1}(i,end) State{2}(i,end)],P,3,'k--')
    end
    
else
    for i = 1:10:length(T)
        P = diag(V(i,1:2)) + V(i,3)*(ones(2)-eye(2));
        
        DrawUniformBounds(E(i,:),P,'k--')
    end
    X = X/length(mus);
end
h1 = plot(X(:,1),X(:,2),'b',lineSpecs{:});
h2 = plot(E(:,1),E(:,2),'r',lineSpecs{:});
h3 = plot(0,0,'k--');


xlabel('x_1',textSpecs{:})
ylabel('x_2',textSpecs{:})
legend([h0, h1, h2, hb, h3],{'f(E[\mu])',['E[f(\mu)] via MC (n=',num2str(length(musMC)),')'],'E[f(\mu)] via PF',boundLabel,'Bound estimate'})

% figure(3)
% surf(State{1},State{2},repmat(mus,n,1),State{3})
% xlabel('x_1',textSpecs{:})
% ylabel('x_2',textSpecs{:})
% zlabel('\mu',textSpecs{:})
%
% figure(4)
% surf(State{1},State{2},State{3},repmat(mus,n,1))
% xlabel('x_1',textSpecs{:})
% ylabel('x_2',textSpecs{:})
% zlabel('p',textSpecs{:})

%
% for i = [1:10:length(T),length(T)]
%     figure(5) hold all plot3(State{1}(i,:), T(i)*ones(1,length(mus)),
%     State{3}(i,:),'*') figure(6) hold all plot3(State{2}(i,:),
%     T(i)*ones(1,length(mus)), State{3}(i,:),'*')
% end

% figure(10000)
% scatter(xf(:,1),xf(:,2),[],musMC)
% hold all
% Pf = diag(V(end,1:2)) + V(end,3)*(ones(2)-eye(2));
% DrawNormalEllipse(E(end,:),Pf,3)

% figure(666) subplot 121 hist(xf(:,1),15) subplot 122 hist(xf(:,2),15)
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


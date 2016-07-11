% Linear Model Predictive Control learning implementation

function u = ModelPredictiveControl()

% Define Test System Model:
n = 4;
m = 2;
p = n;

A = Regularize(rand(n),1);
% B = -2+4*rand(n,m);
B = eye(n);
B(1:2:3,:) = [];
B = B';
% disp('Determining controllable input matrix B')
% tic
% while rank(ctrb(A,B)) < n
%     B = -2+4*rand(n,m);
% end
% disp(['Controllable pair (A,B) determined in ',num2str(toc),' s.'])

C = eye(n);
% C(2:2:n,:) = [];
stri = {'true','false'};
disp(['Observability is ',stri{2-(rank(obsv(A,C)) ==  n)}])

% Set problem info
tf = 25;
x0 = 10*rand(n,1);

% Set weights
Q = eye(n);
R = eye(m);
P = .1*ones(n);

% Set algorithm parameters
dt = 0.1; % How fast the control can change, essentially
N = 10; % Prediction horizon in steps

% Set the parameter structure
pars.n = n;
pars.m = m;
pars.p = p;
pars.N = N;
pars.dt = dt;
pars.tf = tf;
pars.A = A;
pars.B = B;
% pars.K = place(A,B,[1;0.5;0.25;2]);

end

function J = computeCost(x0,U,H,F,Y)
J = U'*H*U + x0'*F*U + x0'*Y*x0;
end

function [G,W,S] = buildConstraints(uLow,uHigh,xLow,xHigh,pars)
I = eye(pars.N);
G = [I;-I];
W = [uHigh*ones(pars.N,1);-uLow*ones(pars.N,1)]; % Control constraints GU <= W

S = 0; % Need to implement state constraints
    
end
    function [H,F,Y] = buildQP(A,B,C,Q,R,P,pars)
    N = pars.N;
    n = pars.n;
    m = pars.m;
    Qb = Tile(Q,pars.N);                            % Create the full matrix
    Qb((end-pars.n):end,(end-pars.n):end) = P;      % Replace the last entry with the final weight matrix
    Rb = Tile(R,pars.N);
    
    S = B;
    T = A;
    Sbar = zeros(N*n,N*m);
    Tbar = zeros(N*n,n);
    for pow = 1:pars.N
        Sbar( 1:n, (1:m*pow) ) = S;
        Tbar((pow-1)*n+(1:n),1:n) = T;
        S = [T*B S];
        T = T*A;
    end
    
    H = (Rb+Sbar'*Qb*Sbar);
    F = 2*Tbar'*Qb*Sbar;
    Y = (Tb'*Qb*Tb)
    
    end

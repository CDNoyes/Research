function srp

clc;
close all

fprintf(['\nA demonstration of the iLQG algorithm '...
    'with SRP dynamics, time as independent variable.\n'...
    'for details see\nTassa, Mansard & Todorov, ICRA 2014\n'...
    '\"Control-Limited Differential Dynamic Programming\"\n'])

% Set full_DDP=true to compute 2nd order derivatives of the
% dynamics. This will make iterations more expensive, but
% final convergence will be much faster (quadratic)

full_DDP = 1;
dist_scale = 3000;
vel_scale = 200;
mass_scale = 7000;

global scale thrust_scale dt
scale = [dist_scale, dist_scale, vel_scale, vel_scale, mass_scale]';
thrust_scale = mass_scale*15;

% set up the optimization problem
DYNCST  = @(x,u,i) entry_dyn_cst(x,u,full_DDP);
tf      = 30;
T       = 100;              % horizon
dt      = tf/T;

x0      = [3000/dist_scale, 3000/dist_scale, -400/vel_scale, -100/vel_scale, 7200/mass_scale]';   % initial state, [DR, Alt, Xdot, Zdot, mass]
u0 = ones(2, T)*7*7200;

Op.lims  = [0 15*7200; 0 15*7200]/thrust_scale;         % wheel angle limits (radians)
Op.plot = 1;               % plot the derivatives as well
Op.maxIter = 50;
Op.parallel = 0;


% === run the optimization!
[x, u, L, Vx, Vxx, cost, trace, stop] = iLQG(DYNCST, x0, u0, Op);

h = x(1,:)*dist_scale/1000;
s = x(2,:)*dist_scale/1000;

vx = x(3,:)*vel_scale;
vz = x(4,:)*vel_scale;
m = x(5,:)*mass_scale;

disp(['hf = ',num2str(h(end)),' km'])
disp(['sf = ',num2str(s(end)),' km'])
disp(['Vx = ',num2str(vx(end)),' m/s'])
disp(['Vz = ',num2str(vz(end)),' m/s'])


figure
plot(s, h)
xlabel('Downrange to go km')
ylabel('Altitude km')
grid on

figure
plot(vx, vz)
% xlabel('Velocity m/s')
% ylabel('Altitude km')
grid on


function y = entry_dynamics(xs, u)
global scale thrust_scale dt

% === states and controls:
% x = [h s vx vz]'
% u = [T]'     = Thrust vector

% constants
% rp = 3396.2e3;

% states
T = u*thrust_scale;

if size(xs,2) == 1
    S = scale;
else
   S = repmat(scale, 1, size(xs, 2)); 
end
x = xs.*S; % need to fix for vectorization
s = x(1,:);
h = x(2,:);
vx = x(3,:);
vz = x(4,:);
m = x(5,:);


g = 3.71;
ve = 295*9.81;

% Derivs
sdot = vx;
hdot = vz;
vxdot = T(1,:)./m;
vzdot = T(2,:)./m - g;
mdot = -sqrt(sum(T.^2, 1))/ve;

xdot = [sdot; hdot; vxdot; vzdot; mdot]./S;            % change in state
y  = xs + xdot.*dt;  % new state


function c = entry_cost(xs, u)
global scale dt thrust_scale

cost_scale = 10000;
% weights


% states
if size(xs,2) == 1
    S = scale;
else
   S = repmat(scale, 1, size(xs, 2)); 
end
T = u*thrust_scale;
x = xs.*S; % need to fix for vectorization
s = x(1,:);
h = x(2,:);
vx = x(3,:);
vz = x(4,:);
m = x(5,:);

g = 3.71;
ve = 295*9.81;

% Derivs
sdot = vx;
hdot = vz;
vxdot = T(1,:)./m;
vzdot = T(2,:)./m - g;
% mdot = -sqrt(sum(T.^2, 1))/ve;

% cost function
% sum of 3 terms:
% lu: quadratic cost on controls
% lf: final cost
% lx: running cost on altitude loss

final = isnan(u(1,:));
u(:,final)  = 0;

% control cost
lu    = 0; %wu*u.^2;

% running cost

% lx = 0.1*sqrt(sum(T/1000.^2, 1)) * dt;
lx = sum(T/100.^2, 1) * dt;
lx = lx + 1e-5 * 2*(s.*sdot + h.*hdot + vz.*vzdot + vx.*vxdot)*dt;


% final cost
if any(final)
    llf      = 0*(0*s.^2 + h.^2) + 1*(vz.^2 + vx.^2); % in real coordinates
    lf       = zeros(1,size(u,2));
    lf(1,final)= llf(1,final);
else
    lf    = 0;
end

% total cost
c     = (lu + lx + lf)/cost_scale;


function [f,c,fx,fu,fxx,fxu,fuu,cx,cu,cxx,cxu,cuu] = entry_dyn_cst(x,u,full_DDP)
% combine dynamics and cost
% use helper function finite_difference() to compute derivatives

if nargout == 2
    f = entry_dynamics(x,u);
    c = entry_cost(x,u);
else
    % state and control indices
    ix = 1:5;
    iu = 6:7;
    
    % dynamics first derivatives
    xu_dyn  = @(xu) entry_dynamics(xu(ix,:),xu(iu,:));
    J       = complex_difference(xu_dyn, [x; u]);
    fx      = J(:,ix,:);
    fu      = J(:,iu,:);
    
    % dynamics second derivatives
    if full_DDP
        N_J = size(J);
        xu_Jcst = @(xu) finite_difference(xu_dyn, xu);
        JJ      = finite_difference(xu_Jcst, [x; u]);
        %         JJ      = reshape(JJ, [4 6 size(J)]); % Original code
        if N_J <= 2
            JJ = reshape(JJ,[ix(end) iu(end) N_J(2)]);
        else
            JJ = reshape(JJ, [ix(end) iu(end) N_J(2) N_J(3)]);
        end
        JJ      = 0.5*(JJ + permute(JJ,[1 3 2 4])); %symmetrize
        fxx     = JJ(:,ix,ix,:);
        fxu     = JJ(:,ix,iu,:);
        fuu     = JJ(:,iu,iu,:);
        
    elseif 0
        N_J = size(J);
        dx = diff([x;u], 1, 2);
        y = diff(J, 1, 3);
        %         I = eye(iu(end));
        B = zeros(N_J(1), N_J(2), N_J(2), N_J(3)); % Hessian
        B(:,:,:,1) = JJ(:,:,:,1);
        for i = 1:(N_J(3)-1) % timesteps
            for j = 1:N_J(1) % each differential equation
                z = y(j,:,i)'-squeeze(B(j,:,:,i))*dx(:,i); % compute (y-Bdx)
                B(j,:,:,i+1) = squeeze(B(j,:,:,i)) + z*z.'/(z.'*dx(:,i));
            end
        end
        E = abs(B-JJ);
        E(isnan(E)) = 0;
        
    else
        [fxx,fxu,fuu] = deal([]);
    end
    
    % cost first derivatives
    xu_cost = @(xu) entry_cost(xu(ix,:),xu(iu,:));
    J       = squeeze(complex_difference(xu_cost, [x; u]));
    cx      = J(ix,:);
    cu      = J(iu,:);
    %     cx_ = gradient(x);
    
    %     cost second derivatives
        xu_Jcst = @(xu) squeeze(complex_difference(xu_cost, xu));    
    
    if 1
        JJ      = finite_difference(xu_Jcst, [x; u]);
        JJ      = 0.5*(JJ + permute(JJ,[2 1 3])); %symmetrize
        cxx     = JJ(ix,ix,:);
        cxu     = JJ(ix,iu,:); % all zeros for alt obj
        cuu     = JJ(iu,iu,:); % all zeros for alt obj

    else   % SR1 estimate
        cxx0     = finite_difference(xu_Jcst, [x(:,1);u(:,1)]);
        JJ      = 0.5*(cxx0 + cxx0.'); %symmetrize
        cxx     = JJ(ix,ix,:);
        N = size(cx);
        dx = diff(x, 1, 2);
        y = diff(cx, 1, 2);
        cxx_qn = zeros(N(1), N(1), N_J(2)); % Hessian
        cxx_qn(:,:,1) = cxx; % initialize with true values
    %     cxx_qn(:,:,1) = eye(N(1)); % initialize with identitity
    
        % how to vectorize this? cannot due to temporal structure 
        for i = 1:(N(2)-1) % timesteps
                z = y(:,i)-cxx_qn(:,:,i)*dx(:,i); % compute (y-Bdx)
                cxx_qn(:,:,i+1) = squeeze(cxx_qn(:,:,i)) + z*z.'/(z.'*dx(:,i));
        end
        cxx = cxx_qn;
        cxu = zeros(N(1), 1, N(2));
        cuu = zeros(1, 1, N(2));
    end
    
    
    [f,c] = deal([]);
end


function J = finite_difference(fun, x, h)
% simple finite-difference derivatives
% assumes the function fun() is vectorized

if nargin < 3
    h = 2^-17;
end

[n, K]  = size(x);
H       = [zeros(n,1) h*eye(n)];
H       = permute(H, [1 3 2]);
X       = pp(x, H);
X       = reshape(X, n, K*(n+1));
Y       = fun(X);
m       = numel(Y)/(K*(n+1));
Y       = reshape(Y, m, K, n+1);
J       = pp(Y(:,:,2:end), -Y(:,:,1)) / h;
J       = permute(J, [1 3 2]);

function J = complex_difference(fun, x)
% simple finite-difference derivatives
% assumes the function fun() is vectorized

h = 0.00001j;

[n, K]  = size(x);
H       = h*eye(n);
H       = permute(H, [1 3 2]);
X       = pp(x, H);
X       = reshape(X, n, K*n);
Y       = fun(X);
m       = numel(Y)/(K*n);
Y       = reshape(Y, m, K, n);
J       = imag(Y)/imag(h);
J       = permute(J, [1 3 2]);

% utility functions: singleton-expanded addition and multiplication
function c = pp(a,b)
c = bsxfun(@plus,a,b);

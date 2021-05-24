% script for creating a database of SRP solutions
clear; clc;
savedir = './data/lower_thrust_weight/';


%% generate samples
N = 4000;
S = lhsdesign(N, 5);
XU = [16000, 1000, 6000, -300,  -50];
XL = [7000,    0,    3000,  -600, -400]; % 800sind(-25)=-340
X0 = floor(XL + S.*(XU-XL));
m0 = 7200;
TW = 4;
T = TW*3.71*m0;
amax = TW*3.71;
% for i = 1:6
%     figure
%     hist((X0(:,i)), 100)
% end


% Remove "bad" states:

% Velocity limit
vmax = 600;
Vmag = sum(X0(:,4:5).^2, 2).^0.5;
k = Vmag < vmax;
disp(sum(k))
X0 = X0(k,:);

% Overshoot limit
% k = overshoot_check(X0(:,1), X0(:,4), m0, T, 290, 0);
% disp(sum(k))
%
% X0 = X0(k,:);

k = overshoot_check(X0(:,3), X0(:,5), m0, T, 290, 3.71);
disp(sum(k))
X0 = X0(k,:);

%%
% X0 = [682, 100, 11375, -785, -77];
% tic;
for i = 1:length(X0) % TODO: try a par for
    
    
    %     x0 = [X0(i, :), m0];
    x0 = [X0(i, 1:4), 0, X0(i,5), m0];
    
  
    if x0(4) < 0
        u = ['n', num2str(-x0(4))];
    else
        u = num2str(x0(4));
    end
    %     if x0(5) < 0
    %         v = ['n', num2str(-x0(5))];
    %     else
    %         v = num2str(x0(5));
    %     end
    if x0(6) < 0
        w = ['n', num2str(-x0(6))];
    else
        w = num2str(x0(6));
    end
    
    
    
    load_guess = 1;
    tf = 0; % for free final time use 0
    
    fname = ['srp_', num2str(x0(1)), '_', num2str(x0(2)), '_', num2str(x0(3)), '_', u, '_', w, '.mat'];
    
    if ~isfile([savedir, fname]) && ~isfile([savedir, 'failed/', fname])  % Verify sample has not been run before
        output = optimize_srp(x0, tf, load_guess);   % Call GPOPS
        sol = output.result.solution.phase(1);
        if output.result.nlpinfo == 1
            save([savedir, fname], 'sol')
            disp(['Solution time = ', num2str(output.totaltime), ' s'])
            disp(['Solution prop = ', num2str(m0+output.result.objective), ' kg'])

            
        else
            save([savedir, 'failed/', fname], 'sol')
            disp(['Solution time = ', num2str(output.totaltime), ' s (did not converge)'])
            
        end
        if ~mod(i, 10)
            clc
            i
        end
    end
end

% timer = toc;
% disp(['Elapsed time for ', num2str(N), ' samples = ', num2str(timer), ' s'])

%%
clear data
data.initial = [];
% data.final = [];
data.tf = [];
data.fuel = [];

files = dir([savedir, '*.mat']);
for i = 1:length(files)
    fname = [savedir, files(i).name];
    load(fname);
    data.initial(end+1,:) = sol.state(1,:);
    %     data.final(end+1,:) = sol.state(end,:);
    data.fuel(end+1) = sol.state(1, end) - sol.state(end, end);
    data.tf(end+1) = sol.time(end);
end
%%
save srp_28k_7200kg.mat -struct data

%% test the iso contours theory
angles = [0, 90, 180, 270];
for i = 1:length(angles)
    x = 20000*cosd(angles(i));
    y = 20000*sind(angles(i));
    v = 500;
    psi = atan2(y, x);
    x0 = [x, y, 5000, -v*cos(psi), -v*sin(psi), -100, m0];
    output = optimize_srp(x0, 0, 0); 
    sol = output.result.solution.phase(1);
    sol.state(1, end) - sol.state(end, end)
end

%%
figure
hist(data.tf)
xlabel('Time of Flight (s)')
figure
hist(data.fuel)
xlabel("Fuel Usage (kg)")

for i = 1:6
    figure
    hist(data.initial(:,i), 500)
    xlabel(num2str(i))
end


% This only makes sense when the downrange was the only variable changing
% [data.initial, Isort] = sortrows(data.initial(:,1));
% data.final = data.final(Isort, :);
% data.fuel = data.fuel(Isort);
% data.tf = data.tf(Isort);
%
% figure
% plot(data.initial(:,1), data.tf)
% xlabel("Initial downrange to target (m)")
% ylabel("Optimal Time of Flight")
% figure
% plot(data.initial(:,1), data.fuel)
% xlabel("Initial downrange to target (m)")
% ylabel("Optimal Fuel Consumed (kg)")

%% Estimate the first order sensitivity of final mass to initial states
x0 = [3200, 0, 3000, -400, 0, -190, 8500];
output = optimize_srp(x0, 0, 0);   % Call GPOPS
nom = output.result.solution.phase.state(end,7);

d = [500, 500, 500, 25, 25, 25]*10;
gplus = zeros(6,1);
gminus = zeros(6,1);
gradient = zeros(6,1);
dm = zeros(6,2);
for i = 1:6
    
    xnew = x0;
    xnew(i) = x0(i) + d(i);
    output = optimize_srp(xnew, 0, 0);   % Call GPOPS
    plus = output.result.solution.phase.state(end,7);
    dm(i,1) = plus-nom;
    
    xnew = x0;
    xnew(i) = x0(i) - d(i);
    output = optimize_srp(xnew, 0, 0);   % Call GPOPS
    minus = output.result.solution.phase.state(end,7);
    dm(i,2) = minus-nom;
    
    gplus(i) = (plus-nom)/d(i);
    gminus(i) = (nom-minus)/d(i);
    gradient(i) = (plus-minus)/(2*d(i));
    
end

% gradient using 500, 25 deltas =
% [-0.098797015592670;-4.528309727902524e-08;-0.185832875842249;
%  -0.237356330662769;-4.059760249219835e-07;-0.677323343770331]

% gradient using 50, 2.5 deltas =
% [-0.104100734339208;-1.216047803609399e-07;-0.187110475373238;
%  -0.088222961046267;9.720388334244490e-08;-0.672907933521492]

% CONCLUSION: Better to use single sided gradient
% Crossrange and side velocity are absolute value looking curves so
% centered difference is roughly zero

%% Redo the above sensitivty for a constant V but changing FPA

Vf = 500;
fpa = -3.5*pi/180;
dfpa = 1 * pi/180;

x0 = @(Vf,fpa) [5500, 0, 2000, -Vf*cos(fpa), 0, Vf*sin(fpa), 8500];

output = optimize_srp(x0(Vf,fpa), 0, 0);   % Call GPOPS
nom = output.result.solution.phase.state(end,7);


dm = zeros(2,2);

xnew = x0(Vf,fpa+dfpa);
output = optimize_srp(xnew, 0, 0);   % Call GPOPS
plus = output.result.solution.phase.state(end,7);
dm(1,1) = plus-nom;

xnew = x0(Vf,fpa-dfpa);
output = optimize_srp(xnew, 0, 0);   % Call GPOPS
minus = output.result.solution.phase.state(end,7);
dm(1,2) = minus-nom;

gplus = (plus-nom)/dfpa;
gminus = (nom-minus)/dfpa;
gradient = (plus-minus)/(2*dfpa);

disp(['gradient is approixmately ',num2str(gradient*pi/180),' kg/deg of FPA'])
%

dv = 20;


xnew = x0(Vf+dv,fpa);
output = optimize_srp(xnew, 0, 0);   % Call GPOPS
plus = output.result.solution.phase.state(end,7);
dm(2,1) = plus-nom;

xnew = x0(Vf-dv,fpa);
output = optimize_srp(xnew, 0, 0);   % Call GPOPS
minus = output.result.solution.phase.state(end,7);
dm(2,2) = minus-nom;


vgplus = (plus-nom)/dv;
vgminus = (nom-minus)/dv;
vgradient = (plus-minus)/(2*dv);

disp(['gradient is approximately ',num2str(vgradient),' kg per m/s of ignition velocity'])

%%
% x0 = [1272.3972 2530.3378 3831.1528  395.2264   -4.4684 -181.7564 8500];
% output = optimize_srp(x0, 0, 0);   % Call GPOPS
% m1 = output.result.solution.phase.state(end,7);
%
% x0 = [1283.  2499.  3759.   397.   -14.5 -170. 8500];
% output = optimize_srp(x0, 0, 0);   % Call GPOPS
% m2 = output.result.solution.phase.state(end,7);

x0 = [ -1.5636342e+03  2.5470960e+03  3.8303374e+03  3.9769870e+02 -1.1210000e+00 -1.8151130e+02 8500];

output = optimize_srp(x0, 0, 0);   % Call GPOPS
m1 = output.result.solution.phase.state(end,7);

x0 = [-1.5700e+03  2.5190e+03  3.7605e+03  3.9800e+02 -3.0000e+00 -1.6800e+02 8500];
output = optimize_srp(x0, 0, 0);   % Call GPOPS
m2 = output.result.solution.phase.state(end,7);

%% Examine the shape of the objective in the vicinity of a good nominal solution

x0 = [ -4000  0  2000  550 0 -100 8500];
dx = [-1000, 1000; 0, 3000; -1000, 1000; -150, 50; 0, 0; -50, 50];
output = optimize_srp(x0, 0, 1);
mf = output.result.solution.phase.state(end,7);
prop0 = 8500-mf;

N = 20;
prop = zeros(6,N);
for i = [1,2,3,4,6]
    xn = x0;
    delta = linspace(dx(i,1), dx(i,2), N);
    for j = 1:N
        xn(i) = x0(i) + delta(j);
        try
            output = optimize_srp(xn, 0, 1);
        catch ME
            output = optimize_srp(xn, 0, 0);
        end
        mf = output.result.solution.phase.state(end,7);
        prop(i,j) = 8500-mf;
    end
    figure(i)
    plot(delta, prop(i,:), '-o')
    hold on
    plot(0, prop0, 'x')
end
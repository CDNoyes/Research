% script for creating a database of SRP solutions

savedir = './data/SRPTable/';

%% generate samples
N = 10000;
S = lhsdesign(N, 6);
XU = [10000, 5000, 4000, 600, 500, 100];
XL = [1000, 0, 1000, 200, -500, -200];
X0 = floor(XL + S.*(XU-XL));

% for i = 1:6
%     figure
%     hist((X0(:,i)), 100)
% end

%%
% x0 = [3200, 0, 3000, -400, 0, -190, 8500];

% tic;
for i = 1:length(X0) % TODO: try a par for 
    x0 = [X0(i, :), 8500];
    
    %     x0 = floor(x0);
    
    if x0(4) < 0
        u = ['n', num2str(-x0(4))];
    else
        u = num2str(x0(4));
    end
    if x0(5) < 0
        v = ['n', num2str(-x0(5))];
    else
        v = num2str(x0(5));
    end
    if x0(6) < 0
        w = ['n', num2str(-x0(6))];
    else
        w = num2str(x0(6));
    end
    
    
    
    load_guess = 0;
    tf = 0; % for free final time use 0
    
    fname = ['srp_', num2str(x0(1)), '_', num2str(x0(2)), '_', num2str(x0(3)), '_', u, '_', v, '_', w, '.mat'];
    
    if ~isfile([savedir, fname])        % Verify sample has not been run before
        output = optimize_srp(x0, 0, 0);   % Call GPOPS
        sol = output.result.solution.phase(1);
        
        if output.result.nlpinfo == 1
            save([savedir, fname], 'sol')
        else
            save([savedir, 'failed/', fname], 'sol')
        end
        if ~mod(i, 50)
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
data.final = [];
data.tf = [];
data.fuel = [];

files = dir([savedir, '*.mat']);
for i = 1:length(files)
    fname = [savedir, files(i).name];
    load(fname);
    data.initial(end+1,:) = sol.state(1,:);
    data.final(end+1,:) = sol.state(end,:);
    data.fuel(end+1) = sol.state(1, end) - sol.state(end, end);
    data.tf(end+1) = sol.time(end);
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
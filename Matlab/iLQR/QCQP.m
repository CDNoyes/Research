function x = QCQP(Q, f, umin, umax)
% Minimize 0.5*x'*Q*x + x'*f  s.t. lower <= ||x|| <= upper
%
%  inputs:
%     Q            - positive definite matrix   (n * n)
%     f            - bias vector                (n)
%     lower        - lower bounds               (1)
%     upper        - upper bounds               (1)

%  outputs:
%     x            - solution                   (n)
%     result       - result type (roughly, higher is better, see below)
%     Hfree        - subspace cholesky factor   (n_free * n_free)
%     free         - set of free dimensions     (n)


c = 0;
% Norm squared constraints umin^2 <= ||u||^2 <= umax^2
H{1} = eye(3); 
H{2} = -eye(3);
k{1} = zeros(3,1);
k{2} = zeros(3,1);
d{1} = -umax^2/2; % negative max^2, times two
d{2} = umin^2/2;  % positive min^2, times two

% On test problems, SQP outperformed IP in time despite using more
% iterations, even with the analytical hessians provided to IP.

% tic;
% options = optimoptions(@fmincon,'Algorithm','interior-point',...
%     'SpecifyObjectiveGradient',true,'SpecifyConstraintGradient',true,...
%     'HessianFcn',@(x,lambda)quadhess(x,lambda,Q,H));
% 
% 
% 
% fun = @(x)quadobj(x,Q,f,c);
% nonlconstr = @(x)quadconstr(x,H,k,d);
% x0 = [0;0;0]; % Column vector
% [x,fval,eflag,output,lambda] = fmincon(fun,x0,...
%     [],[],[],[],[],[],nonlconstr,options);
% disp(x)
% disp(fval)
% disp(output)
% % disp(lambda)
% toc
% tic;
options = optimoptions(@fmincon,'Algorithm','sqp',...
    'SpecifyObjectiveGradient',true,'SpecifyConstraintGradient',true);

fun = @(x)quadobj(x,Q,f,c);
nonlconstr = @(x)quadconstr(x,H,k,d);
x0 = [0;0;0]; % Column vector
[x,fval,eflag,output,lambda] = fmincon(fun,x0,...
    [],[],[],[],[],[],nonlconstr,options);
% disp(x)
% disp(fval)
% disp(output)
% disp(lambda)
% toc

disp('');

function [y,grady] = quadobj(x,Q,f,c)
y = 1/2*x'*Q*x + f'*x + c;
if nargout > 1
    grady = Q*x + f;
end


function [y,yeq,grady,gradyeq] = quadconstr(x,H,k,d)
jj = length(H); % jj is the number of inequality constraints
y = zeros(1,jj);
for i = 1:jj
    y(i) = 1/2*x'*H{i}*x + k{i}'*x + d{i};
end
yeq = [];
    
if nargout > 2
    grady = zeros(length(x),jj);
    for i = 1:jj
        grady(:,i) = H{i}*x + k{i};
    end
end
gradyeq = [];

function hess = quadhess(x,lambda,Q,H)
hess = Q;
jj = length(H); % jj is the number of inequality constraints
for i = 1:jj
    hess = hess + lambda.ineqnonlin(i)*H{i};
end
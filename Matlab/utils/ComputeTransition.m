% COMPUTETRANSITION Computes the state transition matrix of a system along
% a trajectory that has already been designed (or integrated).
%   COMPUTETRANSITION(J,T,X,U,DENSE) Returns the state transition matrix
%   from initial to final states based on the system jacobian J, time
%   vector T and correspond states X, and the control(s) U.
%   DENSE is a flag indicating that all STM matrices from t0 to tf
%   should be output, instead of simply the sensitivity at the final state.



function STM = ComputeTransition(J, time, state, control, dense)
if nargin < 3
    error('Jacobian Function, Independent Variable, and State are required inputs.')
elseif nargin == 3
    control = [];
    dense = 0;
elseif nargin == 4
    dense = 0;
end

[N,n] = size(state);
m = size(control,2); % This allows an empty control to be used to signify no dependence

data.n = n;
data.m = m;
% data.dim = n+m;
data.time = time;
data.state = state;
data.control = control;
data.jacobian = J;

STM_vec0 = reshape(eye(n),[],1);
opt = []; % odeset('AbsTol',1e-9,'RelTol',1e-9);
[~,STM_vec] = ode45(@dynamics, time, STM_vec0,opt,data);

if dense
    % Compute individual STMs here
    STM = cell(1,N);
    for i = 1:N
        STM{i} = reshape(STM_vec(i,:),n,n);
    end
else
    STM = reshape(STM_vec(end,:),n,n);
    
end


end

function dx = dynamics(t,x,data)

% Interp data here
% interp = 'spline';
interp = 'nearest'; % Fast, works fine if your data is on a fine enough grid
state = interp1(data.time,data.state,t,interp);
if data.m
    control = interp1(data.time,data.control,t,interp);
else
    control = [];
end

J = data.jacobian(state,control);
M = reshape(x,data.n,data.n);


% Compute derivative
% dX = J*M;

dx = reshape(J*M,[],1);

end
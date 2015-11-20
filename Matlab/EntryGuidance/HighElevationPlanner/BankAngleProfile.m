%BANKANGLEPROFILE Creates the parameterized profile used in the high
%elevation planner.
%   BANKANGLEPROFILE(T,T1,T2,T3,SIGMA_MIN,SIGMA_MAX) returns the bank angle
%   at time(s) T based on the 3 switching times T1,T2,T3 and the bounds
%   defined by SIGMA_MIN and SIGMA_MAX given in radians.

function Sigma = BankAngleProfile(T,t1,t2,t3,sigma_min, sigma_max)

% Sigma = ones(1,length(T))*-sigma_min;
% Sigma(T>=t1) = sigma_max;
% Sigma(T>=t2) = -sigma_max;
% Sigma(T>=t3) = sigma_min;

%% Slow loop version
% for i = 1:length(T)
%     t = T(i);
%     if t < t1 && t >= 0
%         sigma = -sigma_min;
%     elseif t >= t1 && t < t2
%         sigma = sigma_max;
%     elseif t >= t2 && t < t3
%         sigma = -sigma_max;
%     elseif t >= t3
%         sigma = sigma_min;
%     end
%     Sigma(i) = sigma;
% end


%% Rate-constrained version
dtr = pi/180;
rate = 20*dtr;
% acc = 5 *dtr;

t_inter1 = t1 + (sigma_max+sigma_min)/rate;
t_inter2 = t2 + 2*sigma_max/rate;
t_inter3 = t3 + (sigma_min+sigma_max)/rate;

Sigma = ones(1,length(T))*-sigma_min;
Sigma(T>=t1) =  -sigma_min + rate*(T(T>=t1)-t1);
Sigma(T>=t_inter1) = sigma_max;
Sigma(T>=t2) = sigma_max - rate*(T(T>=t2)-t2);
Sigma(T>=t_inter2) = -sigma_max;
Sigma(T>=t3) = -sigma_max + rate*(T(T>=t3)-t3);
Sigma(T>=t_inter3) = sigma_min;
% if t < t1 && t >= 0
%     sigma = -sigma_min;
% elseif t >= t1 && t <= t_inter1
%     sigma = -sigma_min+rate*(t-t1);
% elseif t > t_inter1 && t < t2
%     sigma = sigma_max;
% elseif t >= t2 && t <= t_inter2
%     sigma = sigma_max - rate*(t-t2);
% elseif t > t_inter2 && t < t3
%     sigma = -sigma_max;
% elseif t > t3 && t <= t_inter3
%     sigma = -sigma_max + rate*(t-t3);
% elseif t > t_inter3
%     sigma = sigma_min;
% end

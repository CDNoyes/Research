function Sigma = BankAngleProfile(T,t1,t2,t3,sigma_min, sigma_max)

Sigma = nan(1,length(T));
for i = 1:length(T)
    t = T(i);
    if t < t1 && t >= 0
        sigma = -sigma_min;
    elseif t >= t1 && t < t2
        sigma = sigma_max;
    elseif t >= t2 && t < t3
        sigma = -sigma_max;
    elseif t >= t3
        sigma = sigma_min;
    end
    Sigma(i) = sigma;
end

% dtr = pi/180;
% rate = 20*dtr;
% acc = 5 *dtr;

% t_inter1 = t1 + (sigma_max+sigma_min)/rate;
% t_inter2 = t2 + 2*sigma_max/rate;
% t_inter3 = t3 + (sigma_min+sigma_max)/rate;
%
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

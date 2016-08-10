function [alpha, k] = computeSMOGains(lambda)
% -lambda is where the three observer poles are placed

% M = [-2,1;-1,0];
% P1 = lyap(M,eye(2)); = [0.5 0.5; 0.5 1.5]
% sqrt(cond(P1)) = 2.4142 so k1 needs to be 2.5 times greater than the
% initial x_1 observation error
C3 = (factorial(3)./(factorial(1:3).*factorial(3-(1:3))));
alpha = C3.*lambda.^(1:3);
k(1) = 5;
k(2:3) = k(1)*[2,1].*lambda.^(1:2);
end
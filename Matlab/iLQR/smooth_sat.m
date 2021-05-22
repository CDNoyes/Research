function y = smooth_sat(x,K)
if nargin == 1
    K = 20;
end
q = 2*x - 1;
y =  0.5/K * log(cosh(K*(q+1))./cosh(K*(q-1))); % saturates between [-1,1]
y = 0.5 + 0.5*y; % map back to [0,1]
end
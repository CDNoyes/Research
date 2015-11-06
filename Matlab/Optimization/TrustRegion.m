function [x,fstar,i,r] = TrustRegion(f,x0)
%Solves a nonlinear, unconstrained minimization problem via trust region
%method. This is meant to be an educational exercise, not an effective
%solver.

x = x0(:);
r = 1; %Initial trust region radius
rmax = 100; %Max radius
imax = 10;
i = 0;
tol = 1e-8;
improvement = 1;
while abs(improvement) > tol && i < imax

fval = f(x);
g = ComplexDiff(f,x);
H = Hessian(f,x);
% D = eye(length(x0)); %Scaling matrix, for problems with parameters of different scale
opt = optimoptions('quadprog','Algorithm','interior-point-convex','Display','off'); %If H is pos def, else maybe use trust region reflective?
opt = optimoptions('quadprog','Algorithm','trust-region-reflective','Display','off');
p = quadprog(H,g,[],[],[],[],-r*ones(size(x)),r*ones(size(x)),[],opt); %Using a square trust region
fstar = f(x+p);
improvement = fstar-fval;
rho = (improvement)/(g*p+0.5*p'*H*p);

if rho >= 0.25
    x = x+p;
    if rho > 0.75
        r = min(2*r,rmax);
    end
else
    r = r/4;
end
i = i+1;
end

end
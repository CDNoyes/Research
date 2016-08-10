function [fIP,jIP,hIP] = InvertedPendulum()

    m = 1; %kg
    l = 0.5; %m
    b = .1;
    I = m*l^2;
    g = 9.81;

    fIP = @(t,x,u) [x(2); (ControlWrapper(t,u) + m*g*l*sin(x(1)) - b*x(2))/I];
    jIP = @(x,u) InvertedPendulumJacobian(x,u,m,l,g,I,b);
    hIP = @(x,u) InvertedPendulumHessian(x,u,m,l,g,I,b);

end

function df = InvertedPendulumJacobian(x,u,m,l,g,I,b)

dfdx = [0, 1; m*g*l*cos(x(1))/I, -b/I];
dfdu = [0;1/I];
df = [dfdx, dfdu];

end


function H = InvertedPendulumHessian(x,u,m,l,g,I,b)

d2fdx2{2} = [-m*g*l*sin(x(1))/I, 0, 0;0,0,0;0,0,0];
d2fdx2{1} = zeros(3);

H = sparse(blkdiag(d2fdx2{:}));
end
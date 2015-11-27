function [fIP,jIP] = InvertedPendulum()

    m = 1; %kg
    l = 0.5; %m
    b = 1;
    I = m*l^2;
    g = 9.81;

    fIP = @(t,x,u) [x(2); (ControlWrapper(t,u) + m*g*l*sin(x(1)) - b*x(2))/I];
    jIP = @InvertedPendulumJacobian;
end

function df = InvertedPendulumJacobian(x,u)

dfdx = [0, 1; m*g*l*cos(x(1)), -b/I];
dfdu = [0;1/I];
df = [dfdx, dfdu];

end
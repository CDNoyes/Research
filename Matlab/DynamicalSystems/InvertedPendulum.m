function [fIP,jIP] = InvertedPendulum()

    fIP = @(t,x,u) [x(2); ControlWrapper(t,u)-4*sin(x(1))];
    jIP = @InvertedPendulumJacobian;
end

function dfdx = InvertedPendulumJacobian(x,u)

dfdx = [0, 1;-4*cos(x(1)) 0];
% dfdu = [0;1];

end
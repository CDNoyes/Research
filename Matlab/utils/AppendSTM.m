% We have to be very particular about the way the inputs are treated. Since
% the function may have state, control, and parameters for which we wish to
% compute sensitivities via STM, the function F should take only one input
% vector X = [x,u,p] and should have trivial dynamics for the constant
% parameters P. Perhaps we can write a wrapper that handles the input
% conversion and trivial dynamics for us. Some thought is required. 

function [fAugmented,xAugmented] = AppendSTM(f,x0,order,jac,hess)

%Order 1 or 2, the order of the state transition matrix stuff to include
nInputs = nargin(f);
switch nInputs
    case 1 % f(t) = f(x) = f(u)
        
    case 2 % f(t,x)
        
    case 3 % f(t,x,u)
        
    case 4 % f(t,x,u,p)
        
    otherwise
        error('Incorrect function input.')
end
xAugmented = [x0(:);reshape(eye(n),[],1)];


end
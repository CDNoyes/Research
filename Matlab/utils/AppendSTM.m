% We have to be very particular about the way the inputs are treated. Since
% the function may have state, control, and parameters for which we wish to
% compute sensitivities via STM, the function F should take only one input
% vector X = [x,u,p] and should have trivial dynamics for the constant
% parameters P. Perhaps we can write a wrapper that handles the input
% conversion and trivial dynamics for us. Some thought is required. 

function [fAugmented,xAugmented] = AppendSTM(f,x0,order,jac,hess)

%Order 1 or 2, the order of the state transition matrix stuff to include
nInputs = nargin(f);


xAugmented = [x0(:);reshape(eye(n),[],1)];


end
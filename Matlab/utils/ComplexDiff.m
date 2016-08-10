%COMPLEXDIFF Complex step differentiation.
%   COMPLEXDIFF(F,X) performs complex differentiation of function F at the
%   point X. Use anonymous functions to create a function with only one
%   input when the original has multiple inputs. E.g.,
%   let g = @(Y) f(Y(1:3),Y(4),Y(5)) when f = @(x,y,z) to compute the
%   gradient of f wrt to x,y,z. If instead, only df/dy (for example) is
%   desired, then use g = @(y) f(x,y,z) where x and z have values in the
%   workspace.
%   
%   NOTE: If the function being differentiated contains the transpose
%   operator, complex differentiation WILL NOT WORK PROPERLY unless the
%   operator .' is used instead. This is because the normal transpose
%   operator returns the complex conjugate for complex numbers.
function dfdx = ComplexDiff(f,x)

n = length(x);
h = 1e-12;

x0 = repmat(x(:),1,n) + h*eye(n)*1j;
for i = 1:n
    dx = f(x0(:,i));
    dfdx(:,i) = imag(dx)/h;
end

end

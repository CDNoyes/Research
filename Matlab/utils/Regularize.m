%REGULARIZE Computes the nearest symmetric, positive (semi)definite
%approximation of A in the sense of the Frobenius norm.
%   REGULARIZE(A,EPSILON) finds a regularized version of the matrix A. When
%   called with only one input argument A, the result X will be positive
%   semidefinite unless A was already positive definite to begin with. When
%   called with a second, positive argument EPSILON, the resulting matrix X
%   is guaranteed to be positive definite with minimum eigenvalue equal to
%   EPSILON. The second output yields the distance between A and X in the
%   sense of Frobenius.

function [X,deltaA,X_inverse] = Regularize(A,epsilon)

if nargin < 2
    epsilon = 0;
else
    assert(epsilon>=0,'Epsilon must be a nonnegative value')
end

B = (A+A')/2; % Symmetric parts of A

%Alternative 1:
% [~,H] = PolarDecomposition(B,1e-12);
% X = (B+H)/2; %Unique positive approximant of A in the Fro norm

%Alternative 2:
[V,D] = eig(B);
d = diag(D);
i_neg = (d < epsilon);
d(i_neg) = epsilon;
X = V*diag(d)*V';
X_inverse = V/diag(d)*V';
if nargout > 1
    d = diag(D);
    C = (A-A')/2; % Skew-symmetric parts of A
    deltaA = sum(d(i_neg).^2)+norm(C,'fro');       
end

end

% function [U,H] = PolarDecomposition(A,tol) %A must be nonsingular
% U = A;
% % n = size(A,1);
% diff = tol+1;
% while diff > tol
%     Ui = inv(U); %can also use U\eye(n) but the speed is comparable
%     gamma = (norm(Ui,1)*norm(Ui,inf)/(norm(U,1)*norm(U,Inf)))^0.25;
%     Unew = 0.5*(gamma*U + Ui'/gamma);
%     diff = norm(Unew-U,1)/norm(Unew,1);
%     U = Unew;
% end
% 
% H = 0.5*(U*A+A*U);
% 
% end
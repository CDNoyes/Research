function xHx = HessianMultiply(H,x,n)

%X' = full(kron(sparse(eye(n)),x'));
%X = repmat(X,n,1)

xHx = full(kron(sparse(eye(n)),x'))*H*repmat(x,n,1);

end
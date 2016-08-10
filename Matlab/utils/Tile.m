%TILE B = tile(A,n)
%Generates a matrix that has the matrix A along its diagonal n times and zeros
%elsewhere. The size of the output matrix is n*size(A). 

function B = tile(A,n)
if n < 2
    error(' n must be greater than 1 ');
end
%Using blkdiag+for loop - Slowest option
% B = A;
% for i = 1:n-1
%     B = blkdiag(B,A);
% end

%Using kron             -Significantly faster
 B = full(kron(sparse(eye(n)),A));
 
%Using blkdiag + more advanced functions - Competitive but not as fast as
%sparse kron
% a = cell(1,n);
% [a{:}] = deal(A);
% B = blkdiag(a{:});

end
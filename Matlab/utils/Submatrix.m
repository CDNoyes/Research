%Implementation that deduces whether the matrix M is a Jacobian or Hessian
%matrix. dim and ind are structures with information about state,control,
%and parameter created using the Dimension function.

function output = Submatrix(M,dim,ind,isScalar)
if ~isstruct(ind) || ~isstruct(dim)
    error('Dim and Ind must be struct inputs, as from Dimension.m')
end
if nargin < 4 || isempty(isScalar)
    isScalar = false;
end
d = size(M);

% isTensor = length(d) == 3; %We never use hessian tensors now
% isHessianCell = isa(M,'cell') ;%Like the second output of Hessian.m
if isScalar
    loopInd = 1;
    isHessianMatrix = d(1) == d(2); %Hessian will be square, Jacobian a vector
else
    loopInd = dim.state:-1:1;
    isHessianMatrix = d(1) ~= dim.state; %Like the first output of Hessian.m
end
output.matrix = M;

if isHessianMatrix
    for i = loopInd
        hxx{i} = M((i-1)*(dim.total)+ind.state,(i-1)*(dim.total)+ind.state);
        hxu{i} = M((i-1)*(dim.total)+ind.state,(i-1)*(dim.total)+ind.control);
        hxp{i} = M((i-1)*(dim.total)+ind.state,(i-1)*(dim.total)+ind.parameter);
        huu{i} = M((i-1)*(dim.total)+ind.control,(i-1)*(dim.total)+ind.control);
        hpp{i} = M((i-1)*(dim.total)+ind.parameter,(i-1)*(dim.total)+ind.parameter);
        hup{i} = M((i-1)*(dim.total)+ind.control,(i-1)*(dim.total)+ind.parameter);
    end
    output.xx = sparse(blkdiag(hxx{:}));
    output.xu = sparse(blkdiag(hxu{:}));
    output.ux = output.xu';
    output.xp = sparse(blkdiag(hxp{:}));
    output.px = output.xp';
    output.uu = sparse(blkdiag(huu{:}));
    output.pp = sparse(blkdiag(hpp{:}));
    output.up = sparse(blkdiag(hup{:}));
    output.pu = output.up';
else %Jacobian
    output.state = M(:,ind.state);
    output.control = M(:,ind.control);
    output.parameter = M(:,ind.parameter);
end


end
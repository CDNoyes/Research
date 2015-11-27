%Implementation that deduces whether the matrix M is a Jacobian or Hessian
%matrix, whether ind is refering to the indices of the matrix or a
%structure for the indices of state,control, and parameter.

function output = Submatrix(M,dim,ind,isScalar)
if ~isstruct(ind)
    error('Ind must be a struct input, as from Dimension.m')
end
if nargin < 4 || isempty(nEquations)
    isScalar = false;
else
    isScalar = true;
end
d = size(M);

% isTensor = length(d) == 3; %We never use hessian tensors now
isHessianCell = isa(M,'cell') ;%Like the second output of Hessian.m
isHessianMatrix = d(1) ~= dim.state; %Like the first output of Hessian.m

output.matrix = M;

if isHessianCell
    %     output.xx
    %     output.xu
    %     output.xp;
    %     output.uu;
    %     output.pp;
    %     output.up;
elseif isHessianMatrix
    for i = dim.state:-1:1
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
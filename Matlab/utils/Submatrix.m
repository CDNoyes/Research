%Implementation that deduces whether the matrix M is a Jacobian or Hessian
%matrix, whether ind is refering to the indices of the matrix or a
%structure for the indices of state,control, and parameter.

function output = Submatrix(M,dim,ind)
if ~isstruct(ind)
    error('Ind must be a struct input, as from Dimension.m')
end

d = size(M);

% isTensor = length(d) == 3; %We never use hessian tensors now
isHessianCell = isa(M,'cell') ;%Like the second output of Hessian.m
isHessianMatrix = d(1) ~= dim.state; %Like the first output of Hessian.m


if isHessianCell
    %     output.xx
    %     output.xu
    %     output.xp;
    %     output.uu;
    %     output.pp;
    %     output.up;
elseif isHessianMatrix
    %H = sparse(blkdiag(h{:}))
    for i = dim.state:-1:1
        hxx{i} = M((i-1)*(dim.total)+ind.state,(i-1)*(dim.total)+ind.state);
    end
    output.xx = sparse(blkdiag(hxx{:}));
    %     output.xu
    %     output.xp;
    %     output.uu;
    %     output.pp;
    %     output.up;
else %Jacobian
    output.state = M(:,ind.state);
    output.control = M(:,ind.control);
    output.parameter = M(:,ind.parameter);
    output.matrix = M;
end


end
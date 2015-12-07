

function sol = ASRE(OCP,A,B,Q,R)

%Handle constant matrices
if ~isa(A,'function_handle') %This should probably never get used
    A = @(x) A;
end
if ~isa(R,'function_handle')
    B = @(x) B;
end
if ~isa(Q,'function_handle')
    Q = @(x) Q;
end
if ~isa(R,'function_handle')
    R = @(x) R;
end

iter = 0;
iterMax = 25;
x0 = OCP.bounds.upper.initialState;
while iter < iterMax
    if ~iter %First iteration is LTI
        K = lqr(A(x0),B(x0),Q(x0),R(x0));
        
    else %Recursive LTV systems
        
        
    end
    
    
    
    iter = iter + 1;
end


end
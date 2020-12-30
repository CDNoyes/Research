function F = LQRGains(Aa, Ba, Qf, Q, R)
%%% A is a vector of nxn matrix, B nxm, R mxm

P = Qf;
F = zeros(size(Ba));

for i = size(Aa,3):-1:2
    
    A = Aa(:,:,i-1);
    B = Ba(:,i-1);
    f = (R + B'*P*B)\(B'*P*A);
    F(:,i-1) = f';
    P = A'*P*A - (A'*P*B)*((R+B'*P*B)\(B'*P*A)) + Q;
    
end
F = -squeeze(F);




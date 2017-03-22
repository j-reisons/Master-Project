function [U,S,V] = rsvd(A,l,q)
%

[m,n] = size(A);

X = randn(n,l);
Y = A*X;
[Q,R] = qr(Y,0);

for j = 1:q
    Y = A'*Q;
    [Q,R] = qr(Y,0);
    Y = A*Q;
    [Q,R] = qr(Y,0);
end

B = Q'*A;
[U,S,V] = svd(B,'econ');
U = Q*U;



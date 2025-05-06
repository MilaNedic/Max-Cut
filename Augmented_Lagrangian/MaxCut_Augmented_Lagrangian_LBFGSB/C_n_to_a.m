function [ A] = C_n_to_a( n, C);
% call [ A] = C_n_to_a( n, C);

m = size( C,1);  % number of sets
hyp=size(C,2);
A = sparse( [],[],[],n*n,m, m*hyp^2);
e = ones(1,hyp);

for i=1:m
   raw = C( i,1:hyp);
   sup1 = abs( raw);
   el1 = sign( raw);
   I = kron(e, sup1);
   J = kron( sup1, e);
   S = kron(el1,el1);
   A0 = sparse( I, J, S, n,n, hyp^2);
   A(:,i) = A0(:);
end
A = A';
function [ A] = C11_to_a( n, C11);
% call [ A] = C5_to_a( n, C5);

m = size( C11,1);  % number of sets
A = sparse( [],[],[],n*n,m, m*121);
e = ones(1,11);

for i=1:m
   raw = C11( i,1:11);
   sup1 = abs( raw);
   el1 = sign( raw);
   I = kron(e, sup1);
   J = kron( sup1, e);
   S = kron(el1,el1);
   A0 = sparse( I, J, S, n,n, 121);
   A(:,i) = A0(:);
end
A = A';

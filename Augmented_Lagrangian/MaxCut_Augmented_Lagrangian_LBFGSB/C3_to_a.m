function [ A] = C3_to_a( n, C3);
% call [ A] = C5_to_a( n, C5);

m = size( C3,1);  % number of sets
A = sparse( [],[],[],n*n,m, m*9);
e = ones(1,3);

for i=1:m
   raw = C3( i,1:3);
   sup1 = abs( raw);
   el1 = sign( raw);
   I = kron(e, sup1);
   J = kron( sup1, e);
   S = kron(el1,el1);
   A0 = sparse( I, J, S, n,n, 9);
   A(:,i) = A0(:);
end
A = A';

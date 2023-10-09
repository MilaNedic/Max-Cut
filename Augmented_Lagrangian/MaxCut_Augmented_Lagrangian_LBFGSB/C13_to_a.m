function [ A] = C13_to_a( n, C13);
% call [ A] = C5_to_a( n, C5);

m = size( C13,1);  % number of sets
A = sparse( [],[],[],n*n,m, m*169);
e = ones(1,13);

for i=1:m
   raw = C13( i,1:13);
   sup1 = abs( raw);
   el1 = sign( raw);
   I = kron(e, sup1);
   J = kron( sup1, e);
   S = kron(el1,el1);
   A0 = sparse( I, J, S, n,n, 169);
   A(:,i) = A0(:);
end
A = A';

function [ A] = C9_to_a( n, C9);
% call [ A] = C5_to_a( n, C5);

m = size( C9,1);  % number of sets
A = sparse( [],[],[],n*n,m, m*81);
e = ones(1,9);

for i=1:m
   raw = C9( i,1:9);
   sup1 = abs( raw);
   el1 = sign( raw);
   I = kron(e, sup1);
   J = kron( sup1, e);
   S = kron(el1,el1);
   A0 = sparse( I, J, S, n,n, 81);
   A(:,i) = A0(:);
end
A = A';

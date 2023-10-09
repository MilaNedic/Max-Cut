function [ A] = C7_to_a( n, C7);
% call [ A] = C5_to_a( n, C5);

m = size( C7,1);  % number of sets
A = sparse( [],[],[],n*n,m, m*49);
e = ones(1,7);

for i=1:m
   raw = C7( i,1:7);
   sup1 = abs( raw);
   el1 = sign( raw);
   I = kron(e, sup1);
   J = kron( sup1, e);
   S = kron(el1,el1);
   A0 = sparse( I, J, S, n,n, 49);
   A(:,i) = A0(:);
end
A = A';

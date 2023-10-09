function [I9, f] = sep_c9( X);
% call: [I5, f] = sep_c9( X);
% use qap

e = ones(9,1);H1 = e*e';
e(1) = -1; H2 = e*e';
e(2) = -1; H3 = e*e';
e(3) = -1; H4 = e*e';
e(4) = -1; H5 = e*e';
trials =1000; %2*size(X,1);%100; % 1000 is good overall % for small graphs this number shouldd not be too big!!
f = zeros( 5*trials, 1); 
I9 = zeros(9, 5*trials); 
for i=1:trials
   [perm,val] = qap_simul2_c( H1, X);
   f(i) = val; 
   I9(:,i) = perm(1:9); 
   
   [perm,val] = qap_simul2_c( H2, X);
   f(trials + i) = val; 
   perm(1) = -perm(1);
   I9(:,trials + i) = perm(1:9); 
   
   [perm,val] = qap_simul2_c( H3, X);
   f(2*trials + i) = val;  
   perm(1) = -perm(1); perm(2) = -perm(2); 
   I9(:,2*trials + i) = perm(1:9); 
   
   [perm,val] = qap_simul2_c( H4, X);
   f(3*trials + i) = val;  
   perm(1) = -perm(1); perm(2) = -perm(2); perm(3) = -perm(3);
   I9(:,3*trials + i) = perm(1:9);
   
   [perm,val] = qap_simul2_c( H5, X);
   f(4*trials + i) = val;  
   perm(1) = -perm(1); perm(2) = -perm(2); perm(3) = -perm(3); perm(4) = -perm(4);
   I9(:,4*trials + i) = perm(1:9);
end
% first check for identical sets
I9keep = I9;
I9 = sort( abs(I9));   % sort columnwise
[dummy, K] = unique(I9', 'rows');
I9 = I9keep(:, K);
f = f(K);
% now remove nonviolations and sort 
K = find(f<0.95); f = f(K); I9=I9(:,K);
[dummy, K] = sort(f);
f = f(K);
I9 = I9(:,K);




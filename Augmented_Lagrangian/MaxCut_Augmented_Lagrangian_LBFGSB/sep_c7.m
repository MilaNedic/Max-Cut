function [I7, f] = sep_c7( X);
% call: [I5, f] = sep_c7( X);
% use qap

e = ones(7,1);H1 = e*e';
e(1) = -1; H2 = e*e';
e(2) = -1; H3 = e*e';
e(3) = -1; H4 = e*e';
trials = 3000;%2*size(X,1);%100; % 1000 is good overall % for small graphs this number shouldd not be too big!!
f = zeros( 4*trials, 1); 
I7 = zeros(7, 4*trials); 
for i=1:trials
   [perm,val] = qap_simul2_c( H1, X);
   f(i) = val; 
   I7(:,i) = perm(1:7); 
   
   [perm,val] = qap_simul2_c( H2, X);
   f(trials + i) = val; 
   perm(1) = -perm(1);
   I7(:,trials + i) = perm(1:7); 
   
   [perm,val] = qap_simul2_c( H3, X);
   f(2*trials + i) = val;  
   perm(1) = -perm(1); perm(2) = -perm(2); 
   I7(:,2*trials + i) = perm(1:7); 
   
   [perm,val] = qap_simul2_c( H4, X);
   f(3*trials + i) = val;  
   perm(1) = -perm(1); perm(2) = -perm(2); perm(3) = -perm(3);
   I7(:,3*trials + i) = perm(1:7);
end
% first check for identical sets
I7keep = I7;
I7 = sort( abs(I7));   % sort columnwise
[dummy, K] = unique(I7', 'rows');
I7 = I7keep(:, K);
f = f(K);
% now remove nonviolations and sort 
K = find(f<0.95); f = f(K); I7=I7(:,K);
[dummy, K] = sort(f);
f = f(K);
I7 = I7(:,K);




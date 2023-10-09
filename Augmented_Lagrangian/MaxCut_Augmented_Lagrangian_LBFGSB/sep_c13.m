function [I11, f] = sep_c13( X);
% call: [I11, f] = sep_c11( X);
% use qap

e = ones(13,1);H1 = e*e';
e(1) = -1; H2 = e*e';
e(2) = -1; H3 = e*e';
e(3) = -1; H4 = e*e';
e(4) = -1; H5 = e*e';
e(5) = -1; H6 = e*e';
e(6) = -1; H7 = e*e';
trials = 110;%2*size(X,1);%100; % 1000 is good overall % for small graphs this number shouldd not be too big!!
f = zeros( 7*trials, 1); 
I11 = zeros(13, 7*trials); 
for i=1:trials
   [perm,val] = qap_simul2_c( H1, X);
   f(i) = val; 
   I11(:,i) = perm(1:13); 
   [perm,val] = qap_simul2_c( H2, X);
   f(trials + i) = val; 
   perm(1) = -perm(1);
   I11(:,trials + i) = perm(1:13); 
   [perm,val] = qap_simul2_c( H3, X);
   f(2*trials + i) = val;  
   perm(1) = -perm(1); perm(2) = -perm(2); 
   I11(:,2*trials + i) = perm(1:13); 
   
   [perm,val] = qap_simul2_c( H4, X);
   f(3*trials + i) = val;  
   perm(1) = -perm(1); perm(2) = -perm(2); perm(3) = -perm(3);
   I11(:,3*trials + i) = perm(1:13);
   
   [perm,val] = qap_simul2_c( H5, X);
   f(4*trials + i) = val;  
   perm(1) = -perm(1); perm(2) = -perm(2); perm(3) = -perm(3); perm(4) = -perm(4);
   I11(:,4*trials + i) = perm(1:13);
   
   [perm,val] = qap_simul2_c( H6, X);
   f(5*trials + i) = val;  
   perm(1) = -perm(1); perm(2) = -perm(2); perm(3) = -perm(3); perm(4) = -perm(4); perm(5) = -perm(5);
   I11(:,5*trials + i) = perm(1:13);
   
   [perm,val] = qap_simul2_c( H7, X);
   f(6*trials + i) = val;  
   perm(1) = -perm(1); perm(2) = -perm(2); perm(3) = -perm(3); perm(4) = -perm(4); perm(5) = -perm(5); perm(6) = -perm(6);
   I11(:,6*trials + i) = perm(1:13);
   
end
% first check for identical sets
I11keep = I11;
I11 = sort( abs(I11));   % sort columnwise
[dummy, K] = unique(I11', 'rows');
I11 = I11keep(:, K);
f = f(K);
% now remove nonviolations and sort 
K = find(f<0.95); f = f(K); I11=I11(:,K);
[dummy, K] = sort(f);
f = f(K);
I11 = I11(:,K);




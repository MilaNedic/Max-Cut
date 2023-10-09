function [I3, f] = sep_c3( X);
% call: [I5, f] = sep_c5( X);
% use qap

e = ones(3,1);H1 = e*e';
e(1) = -1; H2 = e*e';
trials = 100;%2*size(X,1);%100; % 1000 is good overall % for small graphs this number shouldd not be too big!!
f = zeros( 2*trials, 1); 
I3 = zeros(5, 2*trials); 
for i=1:trials
   [perm,val] = qap_simul2_c( H1, X);
   f(i) = val; 
   I3(:,i) = perm(1:5); 
   [perm,val] = qap_simul2_c( H2, X);
   f(trials + i) = val; 
   perm(1) = -perm(1);
   I3(:,trials + i) = perm(1:5); 
end
% first check for identical sets
I3keep = I3;
I3 = sort( abs(I3));   % sort columnwise
[dummy, K] = unique(I3', 'rows');
I3 = I3keep(:, K);
f = f(K);
% now remove nonviolations and sort 
K = find(f<0.95); f = f(K); I3=I3(:,K);
[dummy, K] = sort(f);
f = f(K);
I3 = I3(:,K);




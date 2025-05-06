function [I_n, f] = sep_c_n(X,n,trials);
% n - size of hypermetric ineq.
% trials number of trials for generating hypermetric ineq. to be generated
% use qap
%% call: [I11, f] = sep_c_n(X,n,trials);

if mod(n,2)  %n must ve odd number
    n_2=(n+1)/2;
else
    fprintf('\n\r n must be odd number\n\r')
    return
end
e = ones(n,1);
H=cell(1,n_2);

for ii=1:n_2
    H{ii} = e*e';
    e(ii) = -1;
end;
f = zeros(n_2*trials, 1); 
I_n = zeros(n, n_2*trials); 
for i=1:trials
    for j=1:n_2
        [perm,val] = qap_simul2_c(H{j}, X);
        f(i+(j-1)*trials) = val; 
        perm(1:(j-1)) = -perm(1:(j-1));
        I_n(:,i+(j-1)*trials) = perm(1:n); 
    end
end
% first check for identical sets
I_n_keep = I_n;
I_n = sort( abs(I_n));   % sort columnwise
[dummy, K] = unique(I_n', 'rows');
I_n = I_n_keep(:, K);
f = f(K);
% now remove nonviolations and sort 
K = find(f<0.95); f = f(K); I_n=I_n(:,K);
[dummy, K] = sort(f);
f = f(K);
I_n = I_n(:,K);




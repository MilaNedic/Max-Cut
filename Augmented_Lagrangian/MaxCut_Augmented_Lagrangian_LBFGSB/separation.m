function [B,b] = separation(X,Bold,bold)
% find violated triangle inequalities for primal matrix X
% input:
%       Bold, bold ... old data from previous step
% output:
%       B,b        ... new data, i.e. B(X) <= b
    
n = size(X,1);

% add this number of inequalities
new_ineq = 500; %n*10;

% go to 0-1 model
y = reshape(X, n^2, 1);   % y vector for X

% diagonal elements are not 1 -> need general triangle inequalities
[Tn, g_new] = trianglesep(y, n, new_ineq); 


m = length(g_new);         % number of new ineq.
if m == 0 
    B = Bold;
    b = bold; return    % return with old data
end
if m > 0
    % numbering of cutting planes type is from 1 to 4 instead of 0-3
    Tn(end-m+1:end) = Tn(end-m+1:end) + 1;
end

factor = 1;%sqrt(2/3);

b = [factor*ones(m,1);bold];


[row,column,values] = mc_B(Tn,m,n);

B = sparse(row,column,values,m,n*n,6*m);

B = [factor*B;Bold];


end

function [b, A, X, y,alpha,f] = mcstart(L,alpha,method)
% Setup for max-cut with triangle inequalities
%
% output:
%
%   A,b     ... determines initial set of triangle inequalities
%               A(X) <= b to be included
%   X, y, Z ... starting values for ADMM method obtained from basic
%               SDP relaxation
%
% call:  [b, A, T, X, y, Z] = mcstart(L);


n = size(L,1);
X = eye(n);
y = zeros(n,1);

l  = -inf(n,1);    % there is no lower bound
u  = inf(n,1);      % there is no upper bound
opts    = struct('m', 10, 'printEvery', 0);

opts.factr=0;
%opts.pgtol=0;
opts.pgtol=8e-2;

while true

    switch method
        case 'penalty'
            f = @(y) augLag_penalty(L,y,alpha);
            opts.x0 = y;
            [y,f] = lbfgsb(f,l,u,opts);

            X = 1/alpha * project_W(L - diag(y));

        case 'augLag'
            f = @(y) augLag(L,y,X,alpha);

            opts.x0 = y;
            [y,f] = lbfgsb(f,l,u,opts);

            X = project_W(X + 1/alpha * (L-diag(y)));

            % bound correction
            lambda = min(eig(diag(y) - L));

            f = sum(y) + abs(lambda) * n;

    % update y
    %y0 = y + abs(lambda) * ones(n,1);

    %f = sum(y0);

    end

%f

    % separate triangles
    [A,b] = separation(X,[],[]);

    if size(A,1) == 0
        alpha = alpha * 0.5;
    else
%         if size(A,1) < 50
%             alpha = alpha * 0.8;
%         end
        break;
    end

end

% A = 0;
% b = 0;
% 
% return

% separate triangles
%[A,b] = separation(X,[],[]);

% for i = 1:10
% f = @(y) basic_SDP_dual(L,y,X,alpha);
%     
% opts.x0 = y;
% [y,f] = lbfgsb(f,l,u,opts);
% 
% X = project_W(X + 1/alpha * (L-diag(y)));
% %  
% 
% lambda = min(eig(diag(y) - L));
% 
% % update y
% y0 = y + abs(lambda) * ones(n,1);
% 
% f = sum(y0);




% separate triangles
%[A,b] = separation(X,[],[]);

% if size(A,1) == 0
%     alpha = alpha/2;
% else
%     break;
% end
% 
% end



% print first results
fprintf(' start:   bound: %12.3f\n          triag violation: %.3f,    triag added: %d\n      alpha init: %.3f\n', f, max(A*X(:) - b), length(b), alpha);


end


function [val,grad] = augLag(L,y,X,alpha)
    
    % augmented Lagrangian
    Y = project_W(X + 1/alpha*(L - diag(y)));
    val = sum(y) + alpha/2 * (norm(Y,'fro')^2 - norm(X,'fro')^2);
    
    if nargout > 1
        grad = ones(size(L,1),1) - diag(Y);
    end
    
    

end

function [val,grad] = augLag_penalty(L,y,alpha)
    
    n = size(L,1);
    
    % penalty method 
    Y = project_W(L - diag(y));
    val = sum(y) + 1/(2*alpha) * norm(Y,'fro')^2 + alpha/2 * n^2;
   
    if nargout > 1
        grad = ones(n,1) - 1/alpha*diag(Y);
    end

end
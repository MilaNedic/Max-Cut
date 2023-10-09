function [X, A,b, t, fopt, secs] = main(L,method)
% TODO: MAIN routine for computing SDP relaxation of MAX-CUT strengthen by
% triangle and pentagonal inequalities using Augmented Lagrangian method
% applied to dual SDP
%
% input:    L   ... Laplacian matrix of the graph ( L = diag(Ae)-A )
%           var ... empty (or 3) or 5, determines if triangle or triangle +
%                   pentagonal inequalities are used
%
% output:   X       ... optimal solution of the primal SDP
%           fopt    ... best upper bound
%           secs    ... running time
%
% call: 
%
% for triangles:
%       [X, fopt, secs] = main(L);
% or 
%       [X, fopt, secs] = main(L,3);
%
% for triangles + pentagonal inequalities:
%       [X, fopt, secs] = main(L,5);
%
% ALGORITHM:
% a) compute basic SDP relaxation (interior-point method)
% b) add some number of violated constraints and solve the obtained SDP by
%    ADMM method
% c) purge inactive constraints and go to b)

% TODO: use interior point method for initial relaxtion or use directly
% lbfgs

rng(2020);

% compute first bound and initial set of triangs
tstart = tic;
n = size(L,1);

hyp =3;
alpha = 10;  

[b, A, X, y,alpha,fopt] = mcstart(L,alpha,method);



m = length(b);

z = zeros(m,1); % dual variable to cutting planes       
%y = zeros(n,1); % dual variable to diagonal constraints


maxit = 150;    % number of max iterations
               % (how many times will inequalities be added)
               
             

maxv = 1;      % to initialize stopping condition
maxv5 = 1;

level3 = 0.001;  % violation level for C3 to stop early --> change this to 0.05
level5 = 0.01;  %0.15;


fprintf('        time       bnd         m      purged     maxv3     added3     maxv5     added5     maxv7    added7    maxv9    added9   maxv11   added11\n');

done = 0;

l  = [-inf(n,1); zeros(m,1)];     
u  = inf(n + m,1);          % there is no upper bound


% Request very high accuracy for this test:
opts    = struct('m', 10, 'printEvery', 0);%,'factr', 1e4);
opts.maxIts = 1000;
%opts.factr=1e+4;
%opts.pgtol = 0;

opts.pgtol=8e-2;
scale_tol = 0.9;
min_tol = 1e-5;


% main loop
for cnt = 1:maxit
    
    fold = fopt;
    
    opts.x0 = [y;z];
    m = length(b);
    l  = [-inf(n,1); zeros(m,1)];     
    u  = inf(n + m,1); 
    
    
    switch method
    case 'penalty'
        f = @(var) augLag_penalty(L,var,A,b,alpha);
        
        [var,fopt] = lbfgsb(f,l,u,opts);
        
        y = var(1:n);
        z = var(n+1:end);
        
        
        
        Atz = reshape(A'*z,n,n);
        
        X = 1/alpha * project_W(L - diag(y) - Atz);
       
        
    case 'augLag'
        f = @(var) augLag(L,var,X,A,b,alpha);
        
       
        [var,val] = lbfgsb(f,l,u,opts);
        
        
        y = var(1:n);
        z = var(n+1:end);
        
        Atz = reshape(A'*z,n,n);
        
        X = project_W(X + 1/alpha*(L - diag(y) - Atz));
        %val % lower bound on optimal value of SDP
        
        if alpha > 1
            % bound correction
            %lambda = min(eig(Atz + diag(y) - L))
            
            lambda = max(eig(-(Atz + diag(y) - L)));

            %val - why is this an upper bound??
            fopt = sum(y) + b'*z + lambda * n;
        else
             Z_new = L-diag(y)-Atz;
             [~, y_test, ~, phi] = ipm_mc_pk(Z_new);
             fopt = sum(y) + b'*z + phi;
        end
    end
    
   
    


    % use L-BFGS-B to minimize the dual function
    %f = @(var) dual_function(L,var,X,A,b,alpha);
    
     
%     opts.x0 = [y;z];
%     m = length(b);
%     l  = [-inf(n,1); zeros(m,1)];     
%     u  = inf(n + m,1);          % there is no upper bound
    
%     [var,val,info] = lbfgsb(f,l,u,opts);
%     
%     y = var(1:n);
%     z = var(n+1:end);
%     
%     
%     Atz = reshape(A'*z,n,n);
%     
%      X = project_W(X + 1/alpha*(L - diag(y) - Atz));
    
%     % valid upper bound
%     lambda = min(eig(Atz + diag(y) - L));
%     
%     % update y
%     y0 = y + abs(lambda) * ones(n,1);
%     
%     fopt = sum(y0) + b'*z;
%     
%     fprintf('eig bound: %f\n',fopt);
    
%     Z_new = L-diag(y)-Atz;
% [~, y_test, ~, phi] = ipm_mc_pk(Z_new);
% fopt = sum(y) + b'*z + phi;
    
    %fprintf('sdp bound: %f\n',fopt);
    %

    
    % get primal matrix from dual solution
    


    
    % purge inactive constraints
    old = length(b);
    [b, A, z] = inequality_purge(b, A, z);
    new = length(b);
  
    purged3 = old - new;
    
    m = length(b);       % active constraints

    now = toc(tstart);
    fprintf('%3.0d %8.2f   %10.3f   %4.0d      %4d',cnt, now, fopt, m, purged3);
    
    if ((cnt == maxit) || done);  fprintf('\n'); secs = toc(tstart); return; end
    
    % prune
%     if fopt - lb < 0.99
%         fprintf('\n');
%         break;
%     end
    
    
    % separate new triangles
    old = length(b);            % number of old cutting planes
    [A,b] = separation(X,A,b);
    
    % largest violation of triangles
    maxv = max(A*X(:)-b);   
    fprintf('     %7.3f', maxv);
    
    new_triag = length(b);
    fprintf('       %4d', new_triag - old); % added3
    
    % add pentagonal inequalities if triangles can not reduce bound anymore
    if maxv < 0.3 && (hyp >= 5)  
        old_penta = length(b); % number of old
        [C5,~] = sep_c5(X);
        A5 = C5_to_a(n, C5');
        maxv5 = max(1- A5*X(:));   % largest C5 violation
        
        fprintf('   %7.3f', maxv5);
        
        % new data
        %A = [A; -A5];
        %b = [b;-ones(size(A5,1),1)];
        factor = 1/1;%sqrt(10);
        A = [-factor*A5; A];
        b = [-factor*ones(size(A5,1),1); b];
        
        new_penta = length(b);
        fprintf('       %4d', new_penta - old_penta); % added5
        
    end
    
    if maxv < 0.3 && (hyp >= 7)
        old_hepta = length(b);
        C7 = sep_c7( X);
        A7 = C7_to_a(n, C7');
        maxv7 = max(1- A7*X(:));
        fprintf('   %7.3f', maxv7);   
        
        A = [-A7; A];
        b = [-ones(size(A7,1),1); b];
        
        new_hepta = length(b);
        fprintf('     %4d', new_hepta - old_hepta); % added7

    end
    
    if  maxv < 0.3 && (hyp >= 9)
        old_nona = length(b);
        C9 = sep_c9( X);
        A9 = C9_to_a(n, C9');
        maxv9 = max(1- A9*X(:));
        fprintf('   %7.3f', maxv9);   
        
        A = [-A9; A];
        b = [-ones(size(A9,1),1); b];
        
        new_nona = length(b);
        fprintf('     %4d', new_nona - old_nona); % added9

    end
    
    if (maxv<0.3) && (hyp >= 11)
        old_11 = length(b);
        C11 = sep_c11( X);
        A11 = C11_to_a(n, C11');
        maxv11 = max(1- A11*X(:));
        fprintf('   %7.3f', maxv11);   
        
        A = [-A11; A];
        b = [-ones(size(A11,1),1); b];
        
        new_11 = length(b);
        fprintf('     %4d', new_11 - old_11); % added9

    end
    
    if (maxv<0.3) && (hyp == 13)
        old_13 = length(b);
        C13 = sep_c13( X);
        A13 = C13_to_a(n, C13');
        maxv11 = max(1- A13*X(:));
        fprintf('   %7.3f', maxv11);   
        
        A = [-A13; A];
        b = [-ones(size(A13,1),1); b];
        
        new_13 = length(b);
        fprintf('     %4d', new_13 - old_13); % added9

    end
    
    fprintf('\n');
    
    new = length(b);        % number of new cutting planes
    
    %fprintf('new cuting planes: %d\n',new - old);
    
    % add zero dual multiplier to t for new cuts
    z = [zeros(new-old,1);z]; 
    
    %secs = toc(tstart);
    
    % reduce alpha
    if new_triag - old < 50 %|| (abs(fold - fopt) < 1)
        alpha = alpha * 0.5; % todo: how to decrease?? 0.5 or 0.8
        fprintf('decrease alpha to value %f\n',alpha);
%         hyp = 3; % the  first iter after alpha decreases remove pentagonal cuts
%     else
%         hyp = 5;
         % scale tol
        opts.pgtol = opts.pgtol * scale_tol;
        if opts.pgtol < min_tol
            opts.pgtol = min_tol;
        end
    end 
    
    if (alpha < 1e-6) % && (abs(fold - fopt) < 0.5) % and if bound not close to cut
        break;
    end
    
    % TODO: when decreasing alpha start only with triangles
    % why cant alpha stay the same?
    % look at lbfgs iterations in penalty method
    
end


% 
% for cnt = 1:maxit
% 
        

%     % call ADMM method (change to admm_main_QP for slower quadratic
%     % version
%     [X,y,Z,fopt,t,rho] = admm_main(L,A,b,rho,X,y,Z,t,precision,0);
%       
%     rho = rho/2;
%     
%     precision = precision * 0.8;
%     if (precision < 1e-3)
%         precision = 1e-3;
%     end
%     %gap_precision = gap_precision * 0.7;
%     
%     % purge inactive triangle constraints
%     [b, A, t] = inequality_purge(b, A, t);
%     
%     
%     m = length(b);       % active constraints
% 
%     now = toc(tstart);
%     fprintf('%3.0d %8.2f   %10.6f   %5.0d   %7.3f    %7.3f',cnt, now, fopt, m);
%     
%     
%       
%     if ((cnt == maxit) || done);  fprintf('\n'); secs = toc(tstart); return; end
%     
%     % stopping condition 
%     if var == 3; done = maxv<level3;  end
%     if var == 5; done = maxv5<level5; end
%     
%     % separate new triangles
%     old = length(b);            % number of old cutting planes
%     [A,b] = separation(X,A,b);
%     
%     % largest violation of triangles
%     maxv = max(A*X(:)-b);   
%     fprintf('%7.3f', maxv);
%     
%     % add pentagonal inequalities if triangles can not reduce bound anymore
%     if (maxv<0.2) && (var == 5)      
%         [C5,~] = sep_c5(X);
%         A5 = C5_to_a(n, C5');
%         maxv5 = max(1- A5*X(:));   % largest C5 violation
%         
%         fprintf('   %7.3f', maxv5);
%         
%         % new data
%         %A = [A; -A5];
%         %b = [b;-ones(size(A5,1),1)];
%         factor = 1/sqrt(10);
%         A = [-factor*A5; A];
%         b = [-factor*ones(size(A5,1),1); b];
%         
%     end
%     
%     new = length(b);        % number of new cutting planes
%     
%     %add zero dual multiplier to t for new cuts
%     t = [zeros(new-old,1);t];
%     
%     fprintf('\n');  
%     
%     secs = toc(tstart);
%     
% end

end

function [val,grad] = augLag(L,var,X,A,b,alpha)

    % vector var of varibles. first n correspond to diagonal constraints
    % (free)
    % the other to cutting planes (nonnegative)

    n = size(L,1);
    y = var(1:n);
    z = var(n+1:end);
    
    Atz = reshape(A'*z,n,n);
    
    % inequality auglag
    Y = project_W(X + 1/alpha*(L - diag(y) - Atz));
    val = sum(y) + b'*z + alpha/2 * (norm(Y,'fro')^2 - norm(X,'fro')^2);
    
    
    if nargout > 1
        grad = [ones(n,1) - diag(Y); b - A*Y(:)];
    end
    
    
end

function [val,grad] = augLag_penalty(L,var,A,b,alpha)
    
    n = size(L,1);
    
    y = var(1:n);
    z = var(n+1:end);
    
    Atz = reshape(A'*z,n,n);
    
    % penalty method 
    Y = project_W(L - diag(y) - Atz);
    val = sum(y) + b'*z + 1/(2*alpha) * norm(Y,'fro')^2 + alpha/2 * n^2;
   
    if nargout > 1
        grad = [ones(n,1) - 1/alpha*diag(Y); b - 1/alpha*(A*Y(:))];
    end

end

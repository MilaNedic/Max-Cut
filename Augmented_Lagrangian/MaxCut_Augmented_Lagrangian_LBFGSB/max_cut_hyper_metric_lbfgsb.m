function [X, A,b, fopt,final_output] = max_cut_hyper_metric_lbfgsb(L,method,k,N)
% TODO: MAIN routine for computing SDP relaxation of MAX-CUT strengthen by
% triangle and pentagonal inequalities using Augmented Lagrangian method
% applied to dual SDP
%
% input:    L   ... Laplacian matrix of the graph ( L = diag(Ae)-A )
%           method ... Augm Lang or Augm lagra with penalty term
%           k      ... level of hypermetric inequalties
%           N      ... total number of hypermetric inequalities to be used
% output:   X       ... optimal solution of the primal SDP
%           A,b     ... final hypermetric ineq.
%           fopt    ... best upper bound<
%
%% call: [X, A,b, fopt,final_output] = max_cut_hyper_metric_lbfgsb(L,method,k,N)

rng(2020);

% compute first bound and initial set of triangs
tstart = tic;
n = size(L,1);

hyp =k;
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

fprintf('\n        time       bnd         m      purged');
for i=3:2:k
    if i < 10
        fprintf('     max_v_%1d     added_%1d',i,i);
    elseif i==11 
        fprintf('     max_v_%2d    added_%2d',i,i);
    else
        fprintf('    max_v_%2d    added_%2d',i,i);
    end
end
fprintf('\n')

done = 0;

l  = [-inf(n,1); zeros(m,1)];     
u  = inf(n + m,1);          % there is no upper bound


% Request very high accuracy for this test:
opts    = struct('m', 10, 'printEvery', 0);%,'factr', 1e4);
opts.maxIts = 100; % 100 default
opts.factr=0;
%opts.pgtol = 0;

opts.pgtol=8e-2;
scale_tol = 0.95;
min_tol = 1e-5; %5e-2 default

num_i = ceil(2*N/(hyp-3));

final_output=zeros(ceil((hyp-1)/2),3);

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
    final_output(1,:)=[3 maxv (new_triag - old)];
    
    % add pentagonal inequalities if triangles can not reduce bound anymore
    if maxv < 0.3 
        for hyp_i=5:2:hyp
            old_penta = length(b); % number of old
            [C_n,~] = sep_c_n(X,hyp_i,num_i);
            A_n = C_n_to_a(n, C_n');
            maxv_hyp = max(1- A_n*X(:));   % largest C_n violation
            if isempty(maxv_hyp)
                maxv_hyp = 0;
            end

            fprintf('      %7.3f', maxv_hyp);

            % new data
            %A = [A; -A5];
            %b = [b;-ones(size(A5,1),1)];
            factor = 1/1;%sqrt(10);
            A = [-factor*A_n; A];
            b = [-factor*ones(size(A_n,1),1); b];

            new_penta = length(b);
            fprintf('       %4d', new_penta - old_penta); % added hyp. ineq.
            if floor((hyp_i-1)/2)~=ceil((hyp_i-1)/2) || isempty(hyp_i) || isempty(maxv_hyp) || isempty(new_penta - old_penta)
                keyboard;
            end
            final_output((hyp_i-1)/2,:)=[hyp_i maxv_hyp (new_penta - old_penta)];
        end
    end
    
    fprintf('\n');
    
    new = length(b);        % number of new cutting planes
    
    %fprintf('new cuting planes: %d\n',new - old);
    
    % add zero dual multiplier to z for new cuts
    z = [zeros(new-old,1);z]; 
    
    
    % reduce alpha
    if new_triag - old < 50 || (abs(fold - fopt) < 1)
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
    
    if (alpha < 1e-5) && (abs(fold - fopt) < 0.5) % and if bound not close to 
    % (alpha < 5e-3) && (abs(fold - fopt) < 0.5)
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

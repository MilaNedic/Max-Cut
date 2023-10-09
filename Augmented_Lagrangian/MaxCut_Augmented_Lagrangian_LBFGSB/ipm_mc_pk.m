function [X, y, Z, phi, iter, secs] = ipm_mc_pk(L)
% IPM_MC_PK(primal-dual predictor-corrector interior-point method
% for basic SDP relaxation for Max-Cut) solves:
% 
% max tr(LX), subject to diag(X) = e, X psd
% min r'y, subject to Z = Diag(y) - L psd, y unconstrained
%
% input:  L    ... Laplacian matrix of the graph (1/4 of L for max-cut)
% output: phi  ... optimal value of SDP
%         X    ... optimal primal matrix
%         y    ... optimal dual vector
%         iter ... number of iterations
%         secs ... running time in seconds
%
% call:   [phi, X, y, iter, secs] = ipm_mc(L)

tic;

digits = 6;%5.5;                   % significant digits of phi
n = size(L,1);                % size of problem
b = ones(n,1);                % rhs of linear constraints

% initial positive definite matrices
X = diag(b);
y = 1.1 + sum(abs(L))';
Z = diag(y) - L;


phi = b'*y;                   % initial dual value
psi = L(:)'*X(:);             % initial primal value
mu = Z(:)'*X(:)/(2*n);        % initial complementarity


iter = 0; %iteration count

% while duality gap too large
while (phi - psi) > 1e-2 %max([abs(phi) 1]) * 10^(-digits)
    
    iter = iter + 1;                    % start new iteration
    
    Zi = inv(Z);                        % explicitly compute inv(Z)
    Zi = (Zi + Zi')/2;
      
    
    
    M = Zi.*X;                          % system matrix for dy    
    
    % predictor step (mu = 0) solves: Z * X + Diag(dy1) * X + Z * dX1 = 0
    
    dy1 = M \ ( -b );                   % solve for dy1  
    dX1 = - X - Zi*diag(dy1)*X;         % back substitute for dX 
 
    dX1 = (dX1 + dX1')/2;
    
    
    
    % corrector step solves: diag(dy2)*X + Z*dX2 - mu*I + diag(dy1)*dX1 = 0
    dy2 = M \ (mu*diag(Zi) - (Zi .* dX1)*dy1);
    
    dX2 = mu*Zi - Zi*( diag(dy2) * X + diag(dy1) * dX1);
    
%     if iter == 1
%        disp( - Zi*( diag(dy2) * X + diag(dy1) * dX1)');
%     end
    
    % final steps
    dy = dy1 + dy2;
    dX = dX1 + dX2;
	dX = (dX + dX')/2;          % symmetrise
    dZ = diag(dy);                       % back substitute for dZ
    
    % line search on primal: X  = X + alpha_p * dX  psd matrix
    alpha_p = 1;
    [~,posdef] = chol(X + alpha_p*dX);  % test if positive definite
    
    while posdef > 0
       alpha_p = 0.8 * alpha_p;
       [~,posdef] = chol(X + alpha_p*dX);
    end
    
    if alpha_p < 1
       alpha_p = 0.95 * alpha_p;        % stay away from boundary 
    end
    
    % line search on dual
    alpha_d = 1;
    [~,posdef] = chol(Z + alpha_d*dZ);  
    
    while posdef > 0
       alpha_d = 0.8 * alpha_d;
       [~,posdef] = chol(Z + alpha_d*dZ);
    end
    
    if alpha_d < 1
       alpha_d = 0.95 * alpha_d;        
    end

    % update
    X = X + alpha_p * dX;
    y = y + alpha_d * dy;
    Z = Z + alpha_d * dZ;
    mu = X(:)'*Z(:)/(2*n);
    
    if (alpha_p + alpha_d) > 1.6
        mu = 0.5 * mu;                      % speed up for long steps
    end
    
    if (alpha_p + alpha_d) > 1.9
        mu = 0.2 * mu;                      % speed up for long steps
    end
    
    % objective values
    phi = b'*y;
    psi = L(:)'*X(:);
    
end

secs = toc;


end


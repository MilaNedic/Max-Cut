function [b, A, t]  = inequality_purge(b, A, t)
% Purge inactive triangle and pentagonal constraints. Value of gamma close
% to zero indicates that the corresponding constraint is not binding and we
% can remove it.
%
% call:  [b, A, T, G, gamma]  = bdl_purge( b, A, T, G, gamma);

% determine inactive constraints, given gamma
%maxt = max(t);

% dual variables should be large enough
I = find( t > 1e-5);     
t = t(I);
b = b(I);
A = A(I, :);

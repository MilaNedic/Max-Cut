function [Wp, Wn, V, rank] = project_W( W)
% call: [Wp, Wn, V, rank] = project_W( W);

M = (W + W')/2;      % should be symmetric
n = size( M, 1);

% now compute eigenvalue decomposition of W
[ev, lam] = eig(M); lam = diag(lam);

% compute projection Mp and Mn of M onto PSD and NSD
I = find(lam> 10^-6); j = length(I);
%Computation of Mp: take into account only the "real" positive eigenvalues
evp = zeros( n,j);
for r=1:j;
    ic = I(r); evp(:,r) = ev(:,ic)*sqrt(lam(ic));
end
if j==0; evp = zeros(n,1); end
evpt= evp';
Mp = evp*evpt;
% Computation of V and Mn: take into account only the "real" negative eigenvalues
In = find(lam < -10^-6); jn = length(In);
rank = jn;
pn = zeros(n,jn);
for r=1:jn;
    icn = In(r); pn(:,r) = ev(:,icn)*sqrt(-lam(icn));
end
V = -pn;
pnt = pn';
Mn = -pn*pnt;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Wn = (Mn+Mn')/2; Wp = (Mp+Mp')/2;  % these should be symmetric


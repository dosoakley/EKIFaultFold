function c = CorrMatern(xp,yp,r,nu)
%CorrSpher calculates the correlation matrix for a spherical variogram
%   xp = vector of x locations
%   yp = vector of y locations
%   r = correlation length
%   nu = smoothness parameter
%   c = correlation matrix

N = max(size(xp)); %Note: xp should be only 1 dimensional -- a row or column vector.

%Compute distances. For i==j, dist_mat will be 0, so we don't change it.
dist_mat = zeros(N,N); %Inter-point distances.
for i = 1:N-1
    dists = sqrt((xp(i+1:N)-xp(i)).^2+(yp(i+1:N)-yp(i)).^2);
    dist_mat(i,i+1:N) = dists;
    dist_mat(i+1:N,i) = dists; %Since the distances are the same in either direction.
end

%Compute c
c = (1/(2^(nu-1)*gamma(nu)))*(dist_mat./r).^(nu).*besselk(nu,dist_mat./r);
c(logical(eye(size(c))))=1; %Above equation puts NaN on diagonals, but they should be 1.
end
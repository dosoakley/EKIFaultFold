function r = CorrSpher(xp,yp,Lc)
%CorrSpher calculates the correlation matrix for a spherical variogram
%   xp = vector of x locations
%   yp = vector of y locations
%   Lc = correlation length
%   r = correlation matrix

N = max(size(xp)); %Note: xp should be only 1 dimensional -- a row or column vector.
r = zeros(N,N); %Initialize correlation matrix

%Compute distances. For i==j, dist_mat will be 0, so we don't change it.
dist_mat = zeros(N,N); %Inter-point distances.
for i = 1:N-1
    dists = sqrt((xp(i+1:N)-xp(i)).^2+(yp(i+1:N)-yp(i)).^2);
    dist_mat(i,i+1:N) = dists;
    dist_mat(i+1:N,i) = dists; %Since the distances are the same in either direction.
end

%Compute r
h = dist_mat./Lc;
mask = h<1.0;
r(mask)=1.0+0.5*(-3*h(mask) + h(mask).^3);
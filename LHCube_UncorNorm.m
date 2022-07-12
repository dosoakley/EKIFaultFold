function x = LHCube_UncorNorm(mu,sigma,Ns,varargin)
%Draw Latin hypercube samples normally distributed variables,
%assuming no correlation between the variables.
%This is similar to Matlab's lhsnorm function, but avoids the need for the
%full covariance matrix when working with a purely diagonal covariance
%matrix. Also, the returned x has dimensions nparameters by nsamples, which
%is opposite of lhsnorm but is what I need for this program. And the
%samples are adjusted after sampling to have exactly the desired mean and
%standard deviation.
%Arguments:
%   mu = mean of each variable (1 by Nx element vector).
%   sigma = standard deviation of each variable (1 by Nx element vector).
%   Ns = Number of samples
%   vargin = Name value pairs for lhsdesign.

if iscolumn(mu)
    mu = mu';
end
if iscolumn(sigma)
    sigma = sigma';
end

Nx = length(mu);

if nargin>3
    x = lhsdesign(Ns,Nx,varargin{:}); %Get Latin hypercube samples in the range [0,1].
else
    x = lhsdesign(Ns,Nx);
end

% Transform each column back to the desired marginal distribution,
% maintaining the ranks (and therefore rank correlations) from the
% original random sample
for i=1:Nx
   x(:,i) = norminv(x(:,i),mu(i), sigma(i));
end

%Transpose x.
x = x';

end
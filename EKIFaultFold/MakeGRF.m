function z = MakeGRF(xp,yp,params,zi,variogram,method)
%Create a realization from a Gaussian random field (GRF), using either the 
%Karhunen-Loeve expansion or the Cholesky decomposition.
%Arguments:
%xp,yp = Points at which the GRF should be evaluated.
%params = Variogram parameter. For spherical: range. For Matern: range and smoothness (nu).
%zi = A vector of random numbers drawn from a standard normal distribution.
%variogram = Covariance model to use. Options: 'spherical','matern'
%method = Method for drawing samples. Options: 'kl', 'chol'.
%Notes:
%1) The grid z will be returned in meshgrid format.
%2)The variance is assumed to be 1 here. The result can be multiplied by
%   sqrt(variance) to change that.

switch variogram
    case 'spherical'
        range = params(1);
        c = CorrSpher(xp(:),yp(:),range);
        switch method
            case 'kl'
                [V,D] = eig(c);
                [d,ind] = sort(diag(D),'descend');
                z = V(:,ind)*(sqrt(d).*zi);
            case 'chol'
                L = chol(c,'lower');
                z = L*zi;
            case default
                error('Unrecognized method.')
        end
        z = reshape(z,size(xp));
    case 'matern'
        [range,nu] = deal(params(1),params(2));
        c = CorrMatern(xp(:),yp(:),range,nu);
        switch method
            case 'kl'
                [V,D] = eig(c);
                [d,ind] = sort(diag(D),'descend');
                d(d<0) = 0; %This shouldn't happen, and when it does they tend to be very small, I think due to machine error.
                z = V(:,ind)*(sqrt(d).*zi);
            case 'chol'
                [L,p] = chol(c,'lower');
                if p>=0 %Not positive definite.
                    %c should always be positive definite, but small
                    %machine precision errors for the Matern sometimes mean
                    %it's not.
                    [V,d] = eig(c,'vector');
                    d(d<=1e-10) = 1e-10; %This seems to be about the lowest I can go without chol crashing.
                    c = V*diag(d)/V;
                    L = chol(c,'lower');
                end
                z = L*zi;
            case default
                error('Unrecognized method.')
        end
        z = reshape(z,size(xp));
   
    case default
        error('Unrecognized variogram for gridded GRF.')
end

end
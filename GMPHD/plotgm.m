function [h] = plotgm(mu, sigma, x, y, varargin)
% plotgm(mu, sigma, x, y, varargin) Plot a 2D gaussian mixture
% Inputs:
%   mu          Means of the mixture 2xN, [mu_x, mu_y]
%   sigma       Covariances of the mixture 2x2xN
%   x           X points to plot at 1xK
%   y           Y points to plot at 1xL

# TODO this assumes the gaussian mixture is normalized, which
# isn't necessarily true for an RFS, need to account for weights  
n_gaussians = size(mu, 1);

[X,Y] = meshgrid(x, y);
p = zeros(length(y) * length(x), 1);
for g = 1:n_gaussians
    p = p + mvnpdf([X(:) Y(:)], mu(g, :), sigma(:, :, g));
end

p = reshape(p,size(X));
h = pcolor(X,Y,p);
shading interp;

end

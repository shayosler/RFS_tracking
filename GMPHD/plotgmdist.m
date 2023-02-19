function [h] = plotgmdist(dist, x, y, varargin)
% plotgmdist(dist, x, y, varargin) Plot a 2D gaussian mixture
% Inputs:
%   dist        Distribution. Struct with fields:
%                 
%   x           X points to plot at 1xK
%   y           Y points to plot at 1xL

n_gaussians = size(dist.mu, 1);

[X,Y] = meshgrid(x, y);
p = zeros(length(y) * length(x), 1);
for g = 1:n_gaussians
    if dist.SharedCovariance
        sigma = dist.Sigma;
    else
        sigma = dist.Sigma(:, :, g);
    end
    p = p + dist.ComponentProportion(g) * mvnpdf([X(:) Y(:)], dist.mu(g, :), sigma);
end

p = reshape(p,size(X));
h = pcolor(X,Y,p);
shading interp;

end

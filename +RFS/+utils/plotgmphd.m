function [h] = plotgmphd(dist, x, y, varargin)
% plotgmphd(dist, x, y, varargin) Plot a 2D gaussian mixture probability
% hypothesis density (PHD)
% Inputs:
%   dist    Gaussian mixture representing the PHD   
%   x       X (Northing) points to plot at 1xK
%   y       Y (Easting) points to plot at 1xL

[X,Y] = meshgrid(x, y);
p = zeros(length(y) * length(x), 1);
for g = 1:dist.J
    p = p + dist.w(g) * mvnpdf([X(:) Y(:)], dist.m(:, g)', dist.P(:, :, g));
end

p = reshape(p,size(X));
%h = pcolor(X,Y,p);
h = pcolor(Y, X, p);
shading interp;

end
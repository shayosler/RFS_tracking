function [h] = plotgmphd(dist, x, y, x_ind, y_ind)
% plotgmphd(dist, x, y, varargin) Plot a 2D gaussian mixture probability
% hypothesis density (PHD)
% Inputs:
%   dist    Gaussian mixture representing the PHD   
%   x       X (Northing) points to plot at 1xK
%   y       Y (Easting) points to plot at 1xL
%   x_ind   Index in the state array of the "x" component
%   y_ind   Index in the state array of the "y" component

if nargin < 4 && size(dist.m, 1) == 2
    x_ind = 1;
    y_ind = 2;
end

[X,Y] = meshgrid(x, y);
p = zeros(length(y) * length(x), 1);
for g = 1:dist.J
    p = p + dist.w(g) * mvnpdf([X(:) Y(:)], dist.m([x_ind y_ind], g)', dist.P([x_ind y_ind], [x_ind y_ind], g));
end

p = reshape(p, size(X));
%h = pcolor(X,Y,p);
h = pcolor(Y, X, p);
shading interp;

end
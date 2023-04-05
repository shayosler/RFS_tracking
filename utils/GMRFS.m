classdef GMRFS
    % Gaussian Mixture Random Finite Set
    properties
        J       % Number of components
        m       % Means NxJ
        P       % Covariances NxNxJ
        w       % Weights Jx1
    end

    methods
        function s = GMRFS(mu, sigma, w)
            if (nargin == 0) || ...
                    (isempty(mu) && isempty(sigma) && isempty(w))
                s.J = 0;
                s.m = [];
                s.P = [];
                return;
            end
            s.J = size(mu, 2);
            N = size(mu, 1);
            if size(sigma, 1) ~= N || size(sigma, 2) ~= N || size(sigma, 3) ~= s.J
                error('Invalid sigma');
            end
            if length(w) ~= numel(w) || length(w) ~= s.J
                error('Invalid p')
            end
            w = reshape(w, s.J, 1);
            s.m = mu;
            s.P = sigma;
            s.w = w;

        end

        function s = plus(first, second)
            if first.J == 0
                s = second;
            elseif second.J == 0
                s = first;
            else
                new_mu = [first.m second.m];
                new_sigma = cat(3, first.P, second.P);
                new_w = [first.w; second.w];
                s = GMRFS(new_mu, new_sigma, new_w);
            end
        end

        function s = times(k, rfs)
            if ~isscalar(k)
                error('Only premultiplication by a scalar is supported')
            end
            s = rfs;
            s.w = s.w .* k;
        end

    end
end
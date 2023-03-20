classdef BGRFS
    % Beta-Gaussian Mixture Random Finite Set
    % beta(x) = x^(s-1)(1-x)^(t-1)
    properties
        J       % Number of components, scalar
        m       % Means of gaussian, NxJ
        P       % Covariances of gaussian, NxNxJ
        w       % Weights, Jx1
        s       % Beta distribution s, Jx1
        t       % Beta distribution t, Jx1
    end

    methods
        function d = GMRFS(mu, sigma, w, s, t)
            if (nargin == 0) || ...
                    (isempty(mu) && isempty(sigma) && isempty(w) && isempty(s) && isempty(t))
                d.J = 0;
                d.m = [];
                d.P = [];
                d.s = [];
                d.t = [];
                return;
            end
            d.J = size(mu, 2);
            N = size(mu, 1);
            if size(sigma, 1) ~= N || size(sigma, 2) ~= N || size(sigma, 3) ~= d.J
                error('Invalid sigma');
            end
            if length(w) ~= numel(w) || length(w) ~= d.J
                error('Invalid p')
            end
            if ~isvector(s) || length(s) ~= N
                error('Invalid s')
            end
            if ~isvector(t) || length(t) ~= N
                error('Invalid t')
            end
            w = reshape(w, d.J, 1);
            d.m = mu;
            d.P = sigma;
            d.w = w;
            d.s = s;
            d.t = t;
        end

        function d = plus(first, second)
            % Add two GMRFS objects. The result is a GMRFS object
            % containing the components from both operands
            if first.J == 0
                d = second;
            elseif second.J == 0
                d = first;
            else
                new_mu = [first.m second.m];
                new_sigma = cat(3, first.P, second.P);
                new_w = [first.w; second.w];
                new_s = [first.s; second.s];
                new_t = [first.t; second.t];
                d = RFS.utils.BGRFS(new_mu, new_sigma, new_w, new_s, new_t);
            end
        end

        function s = times(k, rfs)
            % Multiply by a scalar. All weights are multiplied by the
            % specified factor.
            if ~isscalar(k)
                error('Only premultiplication by a scalar is supported')
            end
            s = rfs;
            s.w = s.w .* k;
        end

    end
end

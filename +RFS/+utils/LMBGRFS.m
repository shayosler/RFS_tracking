classdef LMBGRFS
    % Labeled Multi-Bernoulli RFS with Gaussian spatial probability
    properties
        J       % Number of components, scalar        
        w       % Weights of gaussians, Jx1    
        m       % Means of gaussians, NxJ
        P       % Covariances of gaussians, NxNxJ
        r       % Existence probabilities, Jx1  
        l       % Track labels, Jx1
    end

    methods
        function d = LMBGRFS(w, mu, sigma, r, l)
            if (nargin == 0) || ...
                    (isempty(mu) && isempty(sigma) && isempty(w) && isempty(r) && isempty(l))
                d.J = 0;
                d.m = [];
                d.P = [];
                d.r = [];
                d.l = [];             
                return;
            end
            d.J = size(mu, 2);
            N = size(mu, 1);
            if size(sigma, 1) ~= N || size(sigma, 2) ~= N || size(sigma, 3) ~= d.J
                error('Invalid sigma');
            end
            if length(w) ~= numel(w) || length(w) ~= d.J
                error('Invalid w')
            end
            if ~isvector(r) || length(r) ~= d.J
                error('Invalid r')
            end
            if ~isvector(l) || length(l) ~= d.J
                error('Invalid l')
            end
            w = reshape(w, d.J, 1);
            r = reshape(r, d.J, 1);
            l = reshape(l, d.J, 1);
            d.m = mu;
            d.P = sigma;
            d.w = w;
            d.r = r;
            d.l = l;
        end

        function d = plus(first, second)
            % Add two LMBGRFS objects. The result is a LMBGRFS object
            % containing the components from both operands
            if first.J == 0
                d = second;
            elseif second.J == 0
                d = first;
            else
                new_mu = [first.m second.m];
                new_sigma = cat(3, first.P, second.P);
                new_w = [first.w; second.w];
                new_r = [first.r; second.s];
                new_l = [first.l; second.t];
                d = RFS.utils.LMBGRFS(new_w, new_mu, new_sigma, new_r, new_l);
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

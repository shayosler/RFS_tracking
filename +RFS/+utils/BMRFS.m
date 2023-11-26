classdef BMRFS
    % Beta Mixture Random Finite Set
    % beta(x) = x^(s-1)(1-x)^(t-1)
    properties
        J       % Number of components, scalar
        w       % Weights, Jx1
        s       % Beta distribution s, Jx1
        t       % Beta distribution t, Jx1
    end

    methods
        function d = BMRFS(w, s, t)
            if (nargin == 0) || ...
                    (isempty(w) && isempty(s) && isempty(t))
                d.J = 0;
                d.s = [];
                d.t = [];
                return;
            end
            if ~isvector(w)
                error('Invalid w')
            end
            d.J = length(w);
            if ~isvector(s) || length(s) ~= d.J
                error('Invalid s')
            end
            if ~isvector(t) || length(t) ~= d.J
                error('Invalid t')
            end
            w = reshape(w, d.J, 1);
            s = reshape(s, d.J, 1);
            t = reshape(t, d.J, 1);
            d.w = w;
            d.s = s;
            d.t = t;
        end

        function d = plus(first, second)
            % Add two BMRFS objects. The result is a BMRFS object
            % containing the components from both operands
            if first.J == 0
                d = second;
            elseif second.J == 0
                d = first;
            else
                new_w = [first.w; second.w];
                new_s = [first.s; second.s];
                new_t = [first.t; second.t];
                d = RFS.utils.BMRFS(new_w, new_s, new_t);
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

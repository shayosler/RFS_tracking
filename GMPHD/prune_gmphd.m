function [v_out] = prune_gmphd(v, T, U, Jmax)
%[v_out] = prune_gmphd(v, T, U, Jmax) Prune the gaussian mixture components
% representing a random finite set. Pruning is done by removing any 
% components with weight less than T, then merging any components that are
% closer together than U, and finally taking the Jmax highest weighted
% components if there are still more than Jmax
% Inputs
%   v       Gaussian components representing the RFS
%   T       Truncation threshold
%   U       Merging threshold
%   Jmax    Maximum number of components to have after pruning
%
% Outputs
%   v_out   Pruned gaussian mixture

if Jmax <= 0
    warning('Pruning to empty mixture')
    v_out = GMRFS();
    return
end
if v.J < Jmax
    v_out = v;
    return
end

% Set of weights higher than threshold
% TODO: threshold after merging?
I = v.w > T;
idxs = find(I)';

w = [];
m = [];
P = [];
l = 0;
while any(I)
    l = l + 1;
    L = [];
    % Find all components that are close enough to that
    % Largest component to merge
    [~, i_max] = max(v.w(I));
    j = idxs(i_max);
    dist = zeros(sum(I), 1);
    for i = find(I)'
        dist(i) = (v.m(:, i) - v.m(:, j))' * ( v.P(:, :, i) \ (v.m(:, i) - v.m(:, j)));
        if dist(i) < U
            L = [L i];
        end
    end

    % Merge components
    wl = sum(v.w(L));
    ml = (1/wl) .* sum(repmat(v.w(L)', 2, 1) .* v.m(:, L), 2);
    Pl = zeros(size(v.P(:, :, 1)));
    for i = L
        Pl = Pl + v.w(i)*(v.P(:, :, i) + (ml - v.m(:, i))*(ml - v.m(:, i))');
    end
    Pl = Pl ./ wl;
    w = [w wl];
    m = [m ml];
    P = cat(3, P, Pl);

    % Remove merged components
    I(L) = false;
    idxs = find(I)';
end

if numel(w) > Jmax
    [~, imax] = maxk(w, Jmax);
    w = w(imax);
    m = m(:, imax);
    P = P(:, :, imax);
end

v_out = GMRFS(m, P, w);
end
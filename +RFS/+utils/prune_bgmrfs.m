function [v_out] = prune_bgmrfs(v, T, U, Jmax)
%[v_out] = prune_gmphd(v, T, U, Jmax) Prune the gaussian mixture components
% representing a random finite set. Pruning is done by removing any 
% components with weight less than T, then merging any components that are
% closer together than U, and finally taking the Jmax highest weighted
% components if there are still more than Jmax
% Inputs
%   v       Beta-Gaussian components representing the RFS
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
dim = size(v.m, 1);
while any(I)
    l = l + 1;
    % Find all components that are "close" enough to the
    % highest weighted remaining component to merge
    [~, i_max] = max(v.w(I));
    j = idxs(i_max);
    L = []; % Component indices to merge
    dist = zeros(sum(I), 1);
    for i = find(I)'
        dist(i) = (v.m(:, i) - v.m(:, j))' * ( v.P(:, :, i) \ (v.m(:, i) - v.m(:, j)));
        if dist(i) < U
            L = [L i];
        end
    end

    % Merge components
    wl = sum(v.w(L));
    ml = (1/wl) .* sum(repmat(v.w(L)', dim, 1) .* v.m(:, L), 2);
    Pl = zeros(size(v.P(:, :, 1)));
    for i = L
        Pl = Pl + v.w(i)*(v.P(:, :, i) + (ml - v.m(:, i))*(ml - v.m(:, i))');
    end
    Pl = Pl ./ wl;
    w = [w wl];
    m = [m ml];
    P = cat(3, P, Pl);
    sigsq_beta = v.s(L) .* v.t(L) ./((v.s(L) + v.t(L)^2 .* (v.s(L) + v.t(L) + 1)))    ;
    mu_beta = v.s(L) ./ (v.s(L) + v.t(L));
    sigsq_beta_merged = 1 / wl * (v.w(L) * sigsq_beta');
    mu_beta_merged = 1 / wl * (v.w(L) * mu_beta');
    s = 
    t = 

    % Remove merged components
    I(L) = false;
    idxs = find(I)';
end

% If there are still more than Jmax components, remove the lowest weighted
% ones
if numel(w) > Jmax
    [~, imax] = maxk(w, Jmax);
    w = w(imax);
    m = m(:, imax);
    P = P(:, :, imax);
end

% Renormalize weights so that the total weight stays constant
sum_unpruned = sum(v.w);
sum_pruned = sum(w);
w = (sum_unpruned / sum_pruned) .* w;

v_out = RFS.utils.BGMRFS(w, m, P, s, t);
end
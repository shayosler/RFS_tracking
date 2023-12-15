function [v_out] = prune_bgmrfs(v, T, U, Jmax)
%[v_out] = prune_bgmrfs(v, T, U, Jmax) Merge and prune the components of a
% beta-gaussian mixture representing a representing a random finite set. 
% Pruning is done by removing any components with weight less than T, then 
% merging components with a Hellinger distance less than U, and finally 
% taking the Jmax highest weighted components if there are still more than 
% Jmax.
%
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
    v_out = BGMRFS();
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

% Pre-allocate mixture components. Shrink to size later
w = zeros(size(v.w));
m = zeros(size(v.m));
P = zeros(size(v.P));
s = zeros(size(v.s));
t = zeros(size(v.t));
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
    hdist = zeros(sum(I), 1);
    for i = find(I)'
        % FIXME: for some reason this distance calculation 
        % (ie using Hellinger distance) is not giving
        % useful values. Everything has distance very close to 1
        %hdist(i) = RFS.utils.hellinger_bg(v.m(:, i), v.P(:, :, i), v.s(i), v.t(i), v.m(:, j), v.P(:, :, j), v.s(j), v.t(j));
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
    Pl = RFS.utils.make_symmetric(Pl);

    sigsq_beta = v.s(L) .* v.t(L) ./( (v.s(L) + v.t(L)).^2 .* (v.s(L) + v.t(L) + 1) );
    mu_beta = v.s(L) ./ (v.s(L) + v.t(L));
    sigsq_beta_merged = 1 / wl * (v.w(L)' * sigsq_beta);
    mu_beta_merged = 1 / wl * (v.w(L)' * mu_beta);

    com = (mu_beta_merged * (1 - mu_beta_merged) / sigsq_beta_merged) - 1;
    sl = com * mu_beta_merged;
    tl = com * (1 - mu_beta_merged);

    w(l) = wl;
    [~, m(:, l), P(:, :, l)] = RFS.utils.merge_gm(v.w(L), v.m(:, L), v.P(:, :, L));
    [~, s(l), t(l)] = RFS.utils.merge_bm(v.w(L), v.s(L), v.t(L));
    assert(s(l) == sl);
    assert(t(l) == tl);
    assert(all(m(:, l) == ml));
    assert(all(all(P(:, :, l) == Pl)));
    
    % Remove merged components
    I(L) = false;
    idxs = find(I)';
end

% Shrink w, m, P, s, t to the desired size
% If there are more than Jmax components, take the Jmax highest weighted
if l > Jmax
    num = Jmax;
else
    num = l;
end
[~, imax] = maxk(w, num);
w = w(imax);
m = m(:, imax);
P = P(:, :, imax);
s = s(imax, :);
t = t(imax, :);

% Renormalize weights so that the total weight stays constant
sum_unpruned = sum(v.w);
sum_pruned = sum(w);
w = (sum_unpruned / sum_pruned) .* w;

v_out = RFS.utils.BGMRFS(w, m, P, s, t);
end
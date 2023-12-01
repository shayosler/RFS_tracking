function [v_out] = prune_bmrfs(v, T, U, Jmax)
%[v_out] = prune_bmrfs(v, T, U, Jmax) Merge and prune the components of a
% beta mixture representing a representing a random finite set. Pruning is 
% done by removing any components with weight less than T, then merging 
% components with a Hellinger distance less than U, and finally taking the 
% Jmax highest weighted components if there are still more than Jmax.
%
% Inputs
%   v       Beta components representing the RFS
%   T       Truncation threshold
%   U       Merging threshold
%   Jmax    Maximum number of components to have after pruning
%
% Outputs
%   v_out   Pruned beta mixture

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

% Pre-allocate mixture components. Shrink to size later
w = zeros(size(v.w));
s = zeros(size(v.s));
t = zeros(size(v.t));
l = 0;
while any(I)
    l = l + 1;
    % Find all components that are "close" enough to the
    % highest weighted remaining component to merge
    [~, i_max] = max(v.w(I));
    j = idxs(i_max);
    L = []; % Component indices to merge
    dist = zeros(sum(I), 1);
    for i = find(I)'
        dist(i) = RFS.utils.hellinger_beta(v.s(i), v.t(i), v.s(j), v.t(j));
        if dist(i) < U
            L = [L i];
        end
    end

    % Merge components
    [w(l), s(l), t(l)] = RFS.utils.merge_bm(v.w(L), v.s(L), v.t(L));

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
s = s(imax, :);
t = t(imax, :);

% Renormalize weights so that the total weight stays constant
sum_unpruned = sum(v.w);
sum_pruned = sum(w);
w = (sum_unpruned / sum_pruned) .* w;

v_out = RFS.utils.BMRFS(w, s, t);
end
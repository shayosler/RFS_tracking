function [pos] = resize_fig(fig, w, h, preserve_aspect)

if nargin <= 3
    preserve_aspect = false;
end

% Current state
pos = fig.Position;
old_w = pos(3);
old_h = pos(4);
aspect = old_w / old_h;

if w < 0 && h < 0
    error('At least one of w or h must be postive')
elseif w < 0
    % w < 0, choose it to preserve aspect ratio
    new_h = h;
    new_w = aspect * h;
elseif h < 0
    % h < 0, choose it to preserve aspect ratio
    new_w = w;
    new_h = w / aspect;
elseif preserve_aspect
    % Choose new dimensions that preserve aspect
    % ratio and are smaller than w and h
    error('Not implemented');
else
    new_w = w;
    new_h = h;
end

% Update position, preserve size
pos(3) = new_w;
pos(4) = new_h;
fig.Position = pos;

end
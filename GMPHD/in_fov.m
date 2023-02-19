function [in] = in_fov(x, objs, end_range, l_lim, r_lim)
%[in] = in_fov(x, objs, range, l_lim, r_lim) Check if objects are in
%a sonar's field of view
% Inputs:
%   x       vehicle position 3x1, [n; e; psi]
%   objs    Objects to check 2xM [n; e]
%   end_range   Current sonar end range
%   l_lim   Left limit of the sonar field of view relative to current
%           vehicle heading, degrees
%   r_lim   Right limit of the sonar field of view relative to current
%           vehicle heading, degrees
% Outputs:
%   in      1xM logical array that is true for each object in the FOV

% Trivial case
if isempty(objs)
    in = false(0);
    return;
end

n_obj = size(objs, 2);
% Calculate vector from current position to each object
offset = objs - repmat(x(1:2), [1, n_obj]);
r = sqrt(sum(offset.^2, 1));
% relative bearing
b = atan2d(offset(2, :), offset(1, :)) - x(3);

% Determine which targets are within the field of view
% and which ones are actually detected
in = r < end_range & b > l_lim & b < r_lim;

end

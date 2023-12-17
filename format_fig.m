function [] = format_fig(fig, width)

% Set font sizes on all axes
all_ax = findall(fig,'type','axes');
for ax = all_ax
    % Do plot formatting stuff
    set(ax, 'Fontsize', 18)
end

% Resize plot
if width > 0
    %width = 375; % Desired width, pix
    resize_fig(fig, width, -1);
end


end
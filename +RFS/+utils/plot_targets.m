function [h] = plot_targets(targets, varargin)
% plot_targets Plot simulated targets
% 
% Inputs:
%  targets      Vector of RFS.sim.Target_

pos = targets.get_position();
h = plot(pos(2, :), pos(1, :), varargin{:});
end
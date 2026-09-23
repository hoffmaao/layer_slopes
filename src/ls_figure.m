function f = ls_figure(position)
% LS_FIGURE  An invisible, white, light-themed figure for export.
%
%   f = LS_FIGURE([left bottom width height])
%
% From R2025a MATLAB gives new figures the desktop theme, so on a dark
% desktop the axis text comes out light grey and axes backgrounds black -
% nearly invisible in an exported figure with a white background, and a
% black criterion curve on a black axes disappears outright. Every figure
% this repository writes is made here, so it looks the same on every
% release and every desktop. Older releases have no themes and need
% nothing.

f = figure('Visible', 'off');
if isprop(f, 'Theme')
    f.Theme = 'light';
end
set(f, 'Color', 'w', 'Position', position);
end

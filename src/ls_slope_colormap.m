function cmap = ls_slope_colormap(n)
% LS_SLOPE_COLORMAP  Diverging colormap for signed layer slope.
%
%   cmap = LS_SLOPE_COLORMAP(n)
%
% Cyan - blue - white - red - yellow, matching the sense used for reflector
% slope in Holschuh et al. (2017). White sits exactly at zero, so the sign
% of a slope is readable without consulting the colourbar, and the
% saturated cyan/yellow ends make out-of-range values obvious rather than
% letting them blend into the extremes.

if nargin < 1 || isempty(n), n = 256; end

anchors = [ ...
    0.00  0.90 1.00      % cyan
    0.00  0.25 0.90      % blue
    1.00  1.00 1.00      % white
    0.90  0.10 0.10      % red
    1.00  0.90 0.00];    % yellow

t = linspace(0, 1, size(anchors,1));
tq = linspace(0, 1, n);
cmap = [interp1(t, anchors(:,1), tq, 'pchip').', ...
        interp1(t, anchors(:,2), tq, 'pchip').', ...
        interp1(t, anchors(:,3), tq, 'pchip').'];
cmap = min(1, max(0, cmap));
end

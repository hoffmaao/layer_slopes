function c = ls_cice()
% LS_CICE  EM wave speed in solid ice, m/s.
%
%   c = LS_CICE()
%
% c_air / sqrt(3.15) with a relative permittivity of 3.15, giving
% 1.6903e8 m/s. Firn is faster than this, so depths in the top ~100 m
% derived with it are slightly overestimated; that is a systematic scale
% error on depth, not on slope, because a slope is a ratio.
c = 299792458/sqrt(3.15);
end

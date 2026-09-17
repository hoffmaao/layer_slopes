% CICE_IMPORT  Defines cice, the EM wave speed in ice, in the caller's scope.
%
% This is a script, not a function, so that `cice_import` leaves `cice`
% behind in whoever called it - that is how the original code uses it.
%
% The public SlopeExtraction_Radar release calls cice_import from regrid.m
% but does not ship it (it lives in NDH_MatlabTools), so a clean checkout
% cannot regrid radar data. It is included here to make the repo runnable
% standalone.
%
% cair/sqrt(3.15) with a relative permittivity of 3.15 for solid ice gives
% 1.690e8 m/s. RollingRadon.m hard-codes 1.68e8 in one place; the 0.6%
% difference is well inside the uncertainty on the firn correction, but the
% two should not be mixed silently.

cair = 299792458;
cice = cair/sqrt(3.15);

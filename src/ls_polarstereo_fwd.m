function [x, y] = ls_polarstereo_fwd(lat, lon, hemi)
% LS_POLARSTEREO_FWD  Geodetic -> polar stereographic, WGS-84.
%
%   [x,y] = LS_POLARSTEREO_FWD(lat,lon)        Antarctic (EPSG:3031)
%   [x,y] = LS_POLARSTEREO_FWD(lat,lon,hemi)   hemi 0 = south, 1 = north
%
% South uses standard parallel -71 and central meridian 0 (EPSG:3031);
% north uses +70 and -45 (EPSG:3413). Only along-track distance is taken
% from these coordinates, so the projection affects the result only through
% scale - but choosing the wrong hemisphere folds the track back on itself
% and corrupts the distance axis outright.
%
% Snyder, J.P. (1987), Map Projections - A Working Manual, USGS
% Professional Paper 1395, eqns 21-27, 21-32, 21-33.

if nargin < 3 || isempty(hemi), hemi = 0; end

a = 6378137.0;                  % WGS-84 semi-major axis (m)
e = 0.081819190842621;          % WGS-84 first eccentricity

if hemi == 0
    phi_c = -71; lambda_0 = 0;   sgn = -1;
else
    phi_c =  70; lambda_0 = -45; sgn =  1;
end

% Work in the northern formulation, mirror for the south.
phi      = sgn*double(lat)*pi/180;
lambda   = sgn*double(lon)*pi/180;
phi_c    = sgn*phi_c*pi/180;
lambda_0 = sgn*lambda_0*pi/180;

t   = tan(pi/4 - phi/2)   ./ ((1 - e*sin(phi))  ./(1 + e*sin(phi))  ).^(e/2);
t_c = tan(pi/4 - phi_c/2) ./ ((1 - e*sin(phi_c))./(1 + e*sin(phi_c))).^(e/2);
m_c = cos(phi_c)./sqrt(1 - e^2*sin(phi_c).^2);

rho = a*m_c*t/t_c;

x = sgn*  rho.*sin(lambda - lambda_0);
y = sgn*(-rho.*cos(lambda - lambda_0));
end

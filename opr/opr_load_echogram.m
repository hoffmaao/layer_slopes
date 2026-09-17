function D = opr_load_echogram(data_file, opt)
% OPR_LOAD_ECHOGRAM  Read an OPR / CReSIS echogram for slope extraction.
%
%   D = OPR_LOAD_ECHOGRAM(data_file)
%   D = OPR_LOAD_ECHOGRAM(data_file, opt)
%
% opt fields (all optional)
%   .layer_file        explicit CSARP_layer path ('' = derive it)
%   .surface_layer_id  OPR layer id for the surface (default 1)
%   .bed_layer_id      OPR layer id for the bed     (default 2)
%   .verbose           default true
%
% Returns
%   .power        [nt x nx] detected power, linear (|Data|^2 if complex)
%   .data_form    char, how .power was derived from the product
%   .twtt         [nt x 1]  two-way travel time (s)
%   .dist         [1 x nx]  along-track distance (m)
%   .x,.y         [1 x nx]  polar stereographic coordinates (m)
%   .lat,.lon,.elev,.gps_time
%   .surface_twtt [1 x nx]  surface pick (s)
%   .bed_twtt     [1 x nx]  bed pick (s), NaN where unavailable
%   .bed_source   char
%
% THIS FUNCTION NEVER WRITES TO data_file. The original
% RollingRadon_CReSIS.m opened with
%
%       load(filename); clearvars ...; save(filename)
%
% which wrote the function's own arguments (filename, window, plotter,
% movie) back into the shared OPR product. Any echogram that has been
% through that code carries those stowaway variables; the check below
% reports them rather than silently accepting a modified product.
%
% See also ROLLINGRADON_OPR, OPR_LAYER_PATH

if nargin < 2, opt = struct(); end
if ~isfield(opt,'layer_file'),       opt.layer_file = '';  end
if ~isfield(opt,'surface_layer_id'), opt.surface_layer_id = 1; end
if ~isfield(opt,'bed_layer_id'),     opt.bed_layer_id = 2; end
if ~isfield(opt,'verbose'),          opt.verbose = true;   end

data_file = char(data_file);
if exist(data_file,'file') ~= 2
    error('opr_load_echogram:missingFile','Cannot find "%s".', data_file);
end

S = load(data_file);
for req = {'Data','Time','Latitude','Longitude'}
    if ~isfield(S, req{1})
        error('opr_load_echogram:missingVar', ...
            '"%s" has no variable "%s" - is it an OPR echogram?', ...
            data_file, req{1});
    end
end

stowaways = intersect(fieldnames(S), ...
    {'filename','window','plotter','movie','ant_or_green', ...
     'slopes','slope_x','slope_y','slopegrid','slopegrid_x','slopegrid_y'});
if ~isempty(stowaways)
    warning('opr_load_echogram:contaminatedProduct', ...
        ['"%s" carries non-echogram variables (%s). A previous solver run ' ...
         'wrote back into this data product. The echogram itself is still ' ...
         'usable, but this is no longer a pristine OPR file.'], ...
        data_file, strjoin(stowaways', ', '));
end

D = struct();
D.data_file = data_file;

% CSARP_standard / qlook store DETECTED POWER (real). CSARP_standardphase
% and the polarimetric products keep the complex signal, where power is
% |Data|^2. Getting this wrong puts the dB scale out by a factor of two,
% which silently rescales every threshold expressed in dB.
if isreal(S.Data)
    D.power = double(S.Data);
    D.data_form = 'detected power (real)';
else
    D.power = double(abs(S.Data)).^2;
    D.data_form = 'complex signal, converted to |Data|^2';
end
D.twtt  = double(S.Time(:));
D.lat   = double(S.Latitude(:)');
D.lon   = double(S.Longitude(:)');

[nt, nx] = size(D.power);
if numel(D.twtt) ~= nt
    error('opr_load_echogram:shapeMismatch', ...
        'Time has %d samples but Data has %d rows.', numel(D.twtt), nt);
end
if numel(D.lat) ~= nx
    error('opr_load_echogram:shapeMismatch', ...
        'Latitude has %d points but Data has %d columns.', numel(D.lat), nx);
end

% OPR stores seconds. Detect a legacy microsecond record from the total
% record length, not the sample step - a finely sampled microsecond record
% has a small step too, which is what the original heuristic got wrong.
if (D.twtt(end)-D.twtt(1)) > 1e-3
    warning('opr_load_echogram:twttUnits', ...
        'twtt span is %g; assuming microseconds.', D.twtt(end)-D.twtt(1));
    D.twtt = D.twtt*1e-6;
end

if isfield(S,'Elevation'), D.elev = double(S.Elevation(:)'); else, D.elev = nan(1,nx); end
if isfield(S,'GPS_time'),  D.gps_time = double(S.GPS_time(:)'); else, D.gps_time = nan(1,nx); end

% --- geometry -----------------------------------------------------------
hemi = double(median(D.lat) >= 0);   % 0 south (EPSG:3031), 1 north (3413)
[D.x, D.y] = polarstereo_fwd(D.lat, D.lon, hemi);
D.x = D.x(:)'; D.y = D.y(:)';
D.dist = [0 cumsum(hypot(diff(D.x), diff(D.y)))];

% --- surface ------------------------------------------------------------
if isfield(S,'Surface') && any(isfinite(S.Surface))
    D.surface_twtt = double(S.Surface(:)');
    D.surface_source = 'echogram Surface';
else
    D.surface_twtt = zeros(1,nx);
    D.surface_source = 'assumed zero (no Surface in product)';
end

% --- bed ----------------------------------------------------------------
% Ground-based accumulation-radar standard products routinely ship an
% all-NaN Bottom, so fall back to CSARP_layer before giving up.
D.bed_twtt = nan(1,nx);
D.bed_source = 'none';

if isfield(S,'Bottom') && any(isfinite(S.Bottom))
    D.bed_twtt = double(S.Bottom(:)');
    D.bed_source = 'echogram Bottom';
else
    layer_file = opt.layer_file;
    if isempty(layer_file)
        layer_file = opr_layer_path(data_file);
    end
    if ~isempty(layer_file) && exist(layer_file,'file') == 2
        [bed, surf, src] = local_read_layer(layer_file, D.gps_time, opt);
        if any(isfinite(bed))
            D.bed_twtt = bed;
            D.bed_source = src;
        end
        if ~any(isfinite(D.surface_twtt)) && any(isfinite(surf))
            D.surface_twtt = surf;
            D.surface_source = ['CSARP_layer ' layer_file];
        end
    end
end

if opt.verbose
    fprintf('  %s\n', data_file);
    fprintf('    %d samples x %d traces, %.0f m along track\n', nt, nx, D.dist(end));
    fprintf('    data    : %s\n', D.data_form);
    fprintf('    surface : %s\n', D.surface_source);
    fprintf('    bed     : %s (%d/%d finite)\n', ...
        D.bed_source, sum(isfinite(D.bed_twtt)), nx);
end
end

% ------------------------------------------------------------------------
function [bed, surf, src] = local_read_layer(layer_file, gps_time, opt)
bed = nan(size(gps_time));
surf = nan(size(gps_time));
src = '';

L = load(layer_file);
if ~isfield(L,'twtt') || ~isfield(L,'gps_time')
    return
end
if ~any(isfinite(gps_time))
    return
end

nlayer = size(L.twtt,1);
if isfield(L,'id'), ids = L.id(:)'; else, ids = 1:nlayer; end

row = find(ids == opt.bed_layer_id, 1);
if isempty(row)
    [~, row] = max(median(L.twtt, 2, 'omitnan'));   % deepest layer present
end
bed = local_interp(L.gps_time, L.twtt(row,:), gps_time);
src = sprintf('CSARP_layer id=%g (%s)', ids(row), layer_file);

srow = find(ids == opt.surface_layer_id, 1);
if ~isempty(srow)
    surf = local_interp(L.gps_time, L.twtt(srow,:), gps_time);
end
end

% ------------------------------------------------------------------------
function v = local_interp(src_gps, src_val, tgt_gps)
ok = isfinite(src_gps) & isfinite(src_val);
if nnz(ok) < 2
    v = nan(size(tgt_gps));
    return
end
g = src_gps(ok); w = src_val(ok);
[g, iu] = unique(g, 'stable');
v = interp1(g, w(iu), tgt_gps, 'linear', NaN);
end

% DIAG_LAYERING  What does the layering actually look like, and at what
% depth wavelengths does it live? Choose the conditioning from this, not
% from guesswork.

here = fileparts(mfilename('fullpath'));
addpath(fullfile(here,'..','src'));
addpath(fullfile(here,'..','opr'));

data_file = ['/kucresis/scratch/dataproducts/opr_data/accum/' ...
    '2025_Antarctica_Ground2/CSARP_standard_HH/20260109_02/' ...
    'Data_20260109_02_003.mat'];
out_dir = '/kucresis/scratch/hoffmana_sta/layer_slopes/products/diag';
if exist(out_dir,'dir') ~= 7, mkdir(out_dir); end

D = opr_load_echogram(data_file, struct('verbose',false));
cice_import

% Native-resolution, surface-flattened, NO conditioning at all.
G0 = opr_flatten_grid(D, struct('grid_spacing',0.25, 'z_pad_bed',25, ...
    'smooth_len',0, 'detrend_len',0, 'agc_len',0, ...
    'trace_balance',false, 'verbose',true));

fprintf('native dz = %.3f m, grid dz = %.3f m\n', G0.dz_native, G0.grid_spacing);

% ---- depth-wavelength content of the layering --------------------------
% Take a mid-depth band well inside the ice, detrend only the very long
% wavelengths, and look at where the variance lives.
zi = G0.z >= 150 & G0.z <= 450;
band = G0.raw_db(zi,:);
band = band - mean(band,1);
nz = size(band,1);
dz = G0.grid_spacing;

F = abs(fft(band .* hann(nz), [], 1)).^2;
F = mean(F, 2);
half = 2:floor(nz/2);
freq = (half-1)/(nz*dz);              % cycles per metre
wl = 1./freq;                         % metres per cycle
P = F(half);

% cumulative variance from short to long wavelength
[wls, ord] = sort(wl);
Ps = P(ord);
cum = cumsum(Ps)/sum(Ps);
qs = [0.10 0.25 0.50 0.75 0.90];
fprintf('\ndepth-wavelength distribution of layering variance (150-450 m):\n');
for q = qs
    k = find(cum >= q, 1);
    fprintf('   %3.0f%% of variance below %6.1f m wavelength\n', 100*q, wls(k));
end
[~, kpk] = max(Ps);
fprintf('   peak variance at %.1f m wavelength\n', wls(kpk));

% ---- zoom panel at native resolution -----------------------------------
f = figure('Visible','off','Color','w','Position',[50 50 1700 950]);
tl = tiledlayout(f,2,2,'TileSpacing','compact','Padding','compact');

ax = nexttile(tl);
xi = G0.x >= 200 & G0.x <= 700;
zj = G0.z >= 100 & G0.z <= 300;
imagesc(ax, G0.x(xi), G0.z(zj), G0.raw_db(zj,xi));
colormap(ax,gray); set(ax,'YDir','reverse');
clim(ax, prctile(reshape(G0.raw_db(zj,xi),[],1),[5 99]));
ylabel(ax,'depth (m)'); title(ax,'raw dB, 0.25 m grid');

ax = nexttile(tl);
sub = G0.raw_db(zj,xi);
sub = sub - movmean(sub, round(60/dz), 1);
imagesc(ax, G0.x(xi), G0.z(zj), sub);
colormap(ax,gray); set(ax,'YDir','reverse');
clim(ax, prctile(sub(:),[2 98]));
title(ax,'high-pass 60 m');

ax = nexttile(tl);
sub2 = G0.raw_db(zj,xi);
sub2 = sub2 - movmean(sub2, round(15/dz), 1);
imagesc(ax, G0.x(xi), G0.z(zj), sub2);
colormap(ax,gray); set(ax,'YDir','reverse');
clim(ax, prctile(sub2(:),[2 98]));
xlabel(ax,'distance (m)'); ylabel(ax,'depth (m)');
title(ax,'high-pass 15 m (what the last run used)');

ax = nexttile(tl);
loglog(ax, wls, Ps, 'k-'); grid(ax,'on');
xlabel(ax,'depth wavelength (m)'); ylabel(ax,'power');
title(ax,'variance vs depth wavelength');

exportgraphics(f, fullfile(out_dir,'diag_layering.png'), 'Resolution',110);
close(f);
fprintf('\nwrote %s\n', fullfile(out_dir,'diag_layering.png'));

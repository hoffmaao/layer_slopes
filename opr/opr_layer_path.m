function layer_file = opr_layer_path(data_file)
% OPR_LAYER_PATH  Derive the CSARP_layer path matching an OPR echogram path.
%
%   layer_file = OPR_LAYER_PATH(data_file)
%
% OPR lays products out as
%   <season>/CSARP_<product>/<day_seg>/Data_<day_seg>_<frame>.mat
% so the layer product for an echogram is the same path with the product
% directory swapped for CSARP_layer. Image-specific echograms
% (Data_img_01_...) share the segment's single layer file, so the img
% prefix is dropped.
%
% Returns '' when the path does not match that layout.

layer_file = '';
data_file = char(data_file);

[dirpath, name, ext] = fileparts(data_file);
if isempty(dirpath)
    return
end

[season_dir, day_seg] = fileparts(dirpath);
[root, product] = fileparts(season_dir);

if numel(product) < 6 || ~strncmp(product,'CSARP_',6)
    return
end

name = regexprep(name, '^Data_img_\d+_', 'Data_');

layer_file = fullfile(root, 'CSARP_layer', day_seg, [name ext]);
end

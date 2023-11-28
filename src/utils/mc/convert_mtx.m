root_dir = '../../';
data_dir = [root_dir '../sdp_data/'];

src = [data_dir, 'Gset']; % for G set
file_list = dir([src '/*.mat']);
nfile = numel(file_list);
[~, ~] = mkdir(data_dir, 'Gset_mtx');

for i = 1:nfile
    [~, name, ~] = fileparts(file_list(i).name);

    load([src '/' name], 'A', 'n');
    outfile = [data_dir, 'Gset_mtx/', name, '.mtx'];
    fprintf('processing %s ...\n', name);
    mat2sparse(A, outfile);
end

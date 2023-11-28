root_dir = '../../';
data_dir = [root_dir '../sdp_data/'];

subdirs = {'theta', 'theta_new'};

for idS = 1:numel(subdirs)
    setname = subdirs{idS};
    src = [data_dir, setname];
    file_list = dir([src '/*_Alt.mat']);
    nfile = numel(file_list);
    [~, ~] = mkdir(data_dir, [setname, '_mtx']);
    fprintf('Set: %s\n', setname);
    
    for i = 1:nfile
        [~, name, ~] = fileparts(file_list(i).name);
        
        load([src '/' name '.mat'], 'E', 'n', 'p');
        outfile = [data_dir, setname, '_mtx/', name, '.mtx'];
        fprintf('processing %s ...\n', name);
        
        fid = fopen(outfile, 'w');
        % print n, m
        fprintf(fid, '%d %d\n', n, p);
        % print data
        dat = [E'; ones(1, p)];
        fprintf(fid, '%d %d %d\n', dat);
        
        fclose(fid);
    end
end

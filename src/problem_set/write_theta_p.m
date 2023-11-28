file_dir = '.';
root_dir = './..';

data_dir = [root_dir '/sdp_data/'];

result_dir = '../sdpnalplus_results/results/theta/res/';

test_type = 1;
files = thetaprobs_small;
file_len = length(files);
prob_range = 1:file_len;

table_str = ['function ps = thetaprobs_small_p()', newline, '	ps = [ ...', newline];

test_id = 1;

for i = prob_range
    relpath = files{i};
    [~, basename1, basename2] = fileparts(relpath);
    basename = [basename1, basename2];
	load_path = strcat(result_dir,basename,'_res_',mat2str(test_id),'.mat');
	load(load_path);

	p = ceil(0.8*X_rank+0.2*sqrt(2*m));
	cur_str = ['		',mat2str(p),' , ...% ',basename];
	table_str = [table_str, cur_str, newline];
end
table_str = [table_str, '	];', newline];
table_str = [table_str, 'end', newline];

fid = fopen('thetaprobs_small_p.m','w+');
fprintf(fid,'%s',table_str);


table_str = ['function ps = thetaprobs_big_p()', newline, '	ps = [ ...', newline];
files = thetaprobs_big;
file_len = length(files);
prob_range = 1:file_len;

for i = prob_range
    relpath = files{i};
    [~, basename1, basename2] = fileparts(relpath);
    basename = [basename1, basename2];
	load_path = strcat(result_dir,basename,'_res_',mat2str(test_id),'.mat');
	load(load_path);

	p = ceil(0.8*X_rank+0.2*sqrt(2*m));
	cur_str = ['		',mat2str(p),' , ...% ',basename];
	table_str = [table_str, cur_str, newline];
end
table_str = [table_str, '	];', newline];
table_str = [table_str, 'end', newline];

fid = fopen('thetaprobs_big_p.m','w+');
fprintf(fid,'%s',table_str);
function [rankR] = sub_comp_rank(diag_R1, opts)
	if nargin<2; opts = []; end
	if ~isfield(opts,'thres'); opts.thres = 1e-8; end
	if ~isfield(opts,'gap'); opts.gap = 10; end
	if ~isfield(opts,'head'); opts.head = 1; end
	
	thres = opts.thres;
	gap = opts.gap;

	diag_R1 = abs(diag_R1.*(diag_R1>0))+1e-30; % safeguard for division

	rankR = sum(diag_R1>thres*max(diag_R1));
    
    % a negative value of opts.gap means we do not consider gap
    % when estimating rank of R
    if opts.gap > 0
        if opts.head
            jumps = diag_R1(1:rankR-1)./diag_R1(2:rankR);
            inds = find(jumps>=gap);
            if numel(inds)>0
                rankR = inds(end);
            end
        else
            jumps = diag_R1(rankR:end-1)./diag_R1(rankR+1:end);
            inds = find(jumps>=gap);
            if numel(inds)>0
                rankR = rankR-1+inds(1);
            else
                rankR = length(diag_R1);
            end
        end
    end

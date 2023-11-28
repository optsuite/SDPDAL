news = [];%100
news.name = 'news20';
news.raw_path = '../sdp_data/sparse_PCA/raw/20news_w100.mat';
news.multilabel = 0;

madelon = [];%500
madelon.name = 'madelon';
madelon.raw_path = '../sdp_data/sparse_PCA/raw/madelon.txt';
madelon.multilabel = 1;

delicious = [];%500
delicious.name = 'delicious';
delicious.raw_path = '../sdp_data/sparse_PCA/raw/delicious';
delicious.multilabel = 1;

bibtex = [];%1836
bibtex.name = 'delicious';
bibtex.raw_path = '../sdp_data/sparse_PCA/raw/delicious';
bibtex.multilabel = 1;

dna = [];%180
dna.name = 'dna';
dna.raw_path = '../sdp_data/sparse_PCA/raw/dna.scale.txt';
dna.multilabel = 0;

protein = [];%357
protein.name = 'protein';
protein.raw_path = '../sdp_data/sparse_PCA/raw/protein';
protein.multilabel = 0;

mnist = [];%780
mnist.name = 'mnist';
mnist.raw_path = '../sdp_data/sparse_PCA/raw/mnist';
mnist.multilabel = 0;

gisette = [];%5000
gisette.name = 'gisette';
gisette.raw_path = '../sdp_data/sparse_PCA/raw/gisette_scale';
gisette.multilabel = 0;

usps = [];%256
usps.name = 'usps';
usps.raw_path = '../sdp_data/sparse_PCA/raw/usps';
usps.multilabel = 0;

ccancer = [];%2000
ccancer.name = 'colon_cancer';
ccancer.raw_path = '../sdp_data/sparse_PCA/raw/colon-cancer';
ccancer.multilabel = 0;

datasets = [news,madelon,delicious,bibtex,dna,protein,mnist,ccancer,gisette,usps];

for idx = 1:length(datasets)
    dataset_cur = datasets(idx);
    save_path = sprintf('../sdp_data/sparse_PCA/%s.mat',dataset_cur.name);
    fprintf('%s\n',save_path)
    if strcmp(dataset_cur.name,'news20')
        load(dataset_cur.raw_path)
        data = documents';
        data_mean = mean(data,1);
        data = data-data_mean;
        Sigma = data'*data;
        Sigma = full(Sigma);
    elseif dataset_cur.multilabel
        [b, A ] = libsvmread_ml(dataset_cur.raw_path);
        A_mean = mean(A,1);
        A = A-A_mean;
        Sigma = A'*A;
        Sigma = full(Sigma);
    else
        [b, A ] = libsvmread(dataset_cur.raw_path);
        A_mean = mean(A,1);
        A = A-A_mean;
        Sigma = A'*A;
        Sigma = full(Sigma);
    end
    save(save_path, 'Sigma')
    
end
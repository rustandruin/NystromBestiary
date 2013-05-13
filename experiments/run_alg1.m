% Computes Nystrom extensions using approximate leverage score sampling, where the 
% leverage scores are approximated using Algorithm 1 from 
% "Fast approximation of matrix coherence and statistical leverage"
% on two datasets with appropriately small aspect ratio (tall and thin)
% and stores the output

tallthinmethod = {'tallthinlevscoreapprox'};
datadirname = 'outputs'; % where to put the generated data

p = .9;
q = 6;
numpts = 12;
chunk = 10;
vareps = 1/3;
tol = .01; %1/3
tallthineps = 1;

genlvals = @(startl, endl) ceil(linspace(startl, endl, numpts));


%% Protein

load 'protein_data';
in.A = normalize_kernel_data(X);
clear X;
in.linearkernelflag = 1;
in.k = 10; % captures 42% of the frob norm, doubling to 20 captures only 51%
in.chunk = chunk;
in.vareps = vareps;
in.tol = tol;
in.lvals = genlvals(in.k, 110);
in.p = p;
in.q = q;
in.tallthineps = 1; 
in.tallthinmethod = 'svd'; % Don't forget to change approx_tallthin_levscores to use SVD since Protein is rank-deficient, so QR gives warning

in.descr = '(whitened) Protein linear kernel';
in.datasetbasename = 'tallthinProteinrank10';
in.datasetdir = datadirname;
in.methods = tallthinmethod;

generate_dataset(in);
clear in;

%% SNPS 
load 'snpsdata'
A = A(1:(end-3), 1:size(A,2)-3); % The last three columns and rows are all 0s
in.A = normalize_kernel_data(A);
clear A;
in.linearkernelflag = 1;
in.k = 5; % A is clearly low rank
in.chunk = chunk;
in.vareps = vareps;
in.tol = tol;
in.lvals = genlvals(in.k, 46);
in.p = p;
in.q = q;
in.tallthineps = 1; 
in.tallthinmethod = 'qr';

in.descr = '(whitened) SNPS linear kernel';
in.datasetbasename = 'tallthinSNPSrank5';
in.datasetdir = datadirname;
in.methods = tallthinmethod;

generate_dataset(in);
clear in;


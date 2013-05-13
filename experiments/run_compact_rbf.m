%% Runs Nystrom extensions on the compact RBF datasets and stores 
% results

% specify which methods to use
allmethods = {'simple', 'srft', 'gaussian', 'levscore', 'froblev', ...
    'speclev', 'approxlev'};
% specify the directory to output the results into
datadirname = 'outputs';

p = .9; % for the mixedprobs method
q = 6; % number of experiments to run
numpts = 12; % number of different column samples to take
chunk = 10; % how often to reorthogonalize power method
vareps = 1/3; % ignored, is the epsilon parameter in some methods
tol = .01; % the tolerance to seek the leverage scores to in the power method
cutoffmultiplier = 3; % parameter for the generation of compact rbf matrices

% generates the set of numbers of column samples to be taken
genlvals = @(startl, endl) ceil(linspace(startl, endl, numpts));


%% Abalone, sigma = 1

in.sigma = 1;
in.cutoff = cutoffmultiplier*in.sigma;
in.d = 8; % number of features per observation

load 'abalone_distance_matrix'; 
in.A = generate_compact_RBF_kernel(X, in.sigma, in.d, in.cutoff);
in.sparsity = nnz(in.A)/prod(size(in.A));
clear X;

in.linearkernelflag = 0; 
in.k = 20; % target tank
in.chunk = chunk;
in.vareps = vareps;
in.tol = tol;
in.lvals = genlvals(in.k, 150);
in.p = p;
in.q = q;

in.descr = 'Abalone compact RBF kernel with sigma = 1';
in.datasetbasename = 'Abalonecompactsigma1';
in.datasetdir = datadirname;
in.methods = allmethods;

generate_dataset(in);
clear in;

%% Abalone, sigma = .15

in.sigma = .15; % when sigma = .15, k=20 captures < 18% of variance
in.cutoff = cutoffmultiplier*in.sigma;
in.d = 8;

load 'abalone_distance_matrix'; 
in.A = generate_compact_RBF_kernel(X, in.sigma, in.d, in.cutoff);
in.sparsity = nnz(in.A)/prod(size(in.A));
clear X;

in.linearkernelflag = 0;
in.k = 20;
in.chunk = chunk;
in.vareps = vareps;
in.tol = tol;
in.lvals = genlvals(in.k, 150);
in.p = p;
in.q = q;

in.descr = 'Abalone compact RBF kernel with sigma = .15';
in.datasetbasename = 'Abalonecompactsigmapt15';
in.datasetdir = datadirname;
in.methods = allmethods;

generate_dataset(in);
clear in;

%% Wine, sigma = 1

load 'winequality_distance_matrix';
in.sigma = 1;
in.d = 12;
in.cutoff = cutoffmultiplier*in.sigma;

in.A = generate_compact_RBF_kernel(X, in.sigma, in.d, in.cutoff);
in.sparsity = nnz(in.A)/prod(size(in.A));
clear X;

in.linearkernelflag = 0;
in.k = 20;
in.chunk = chunk;
in.vareps = vareps;
in.tol = tol;
in.lvals = genlvals(in.k, 150);
in.p = p;
in.q = q;

in.descr = 'Wine compact RBF kernel with sigma = 1';
in.datasetbasename = 'Winecompactsigma1';
in.datasetdir = datadirname;
in.methods = allmethods;

generate_dataset(in);
clear in;

%% Wine, sigma = 2.1

load 'winequality_distance_matrix';
in.sigma = 2.1;
in.d = 12;
in.cutoff = cutoffmultiplier*in.sigma;

in.A = generate_compact_RBF_kernel(X, in.sigma, in.d, in.cutoff);
in.sparsity = nnz(in.A)/prod(size(in.A));
clear X;

in.linearkernelflag = 0;
in.k = 20;
in.chunk = chunk;
in.vareps = vareps;
in.tol = tol;
in.lvals = genlvals(in.k, 150);
in.p = p;
in.q = q;

in.descr = 'Wine compact RBF kernel with sigma = 2.1';
in.datasetbasename = 'Winecompactsigma2pt1';
in.datasetdir = datadirname;
in.methods = allmethods;

generate_dataset(in);
clear in;

%% Kin8nm

load 'kin8nm_distance_matrix'
in.sigma = 2.2;
in.d = 9;
in.cutoff = cutoffmultiplier*in.sigma;

in.A = generate_compact_RBF_kernel(kin8nm_dist, in.sigma, in.d, in.cutoff);
in.sparsity = nnz(in.A)/prod(size(in.A));
clear kin8nm_dist;

in.linearkernelflag = 0;
in.k = 40;
in.chunk = chunk;
in.vareps = vareps;
in.tol = tol;
in.lvals = genlvals(in.k, 150);
in.p = p;
in.q = q;

in.descr = 'Kin8nm RBF kernel with sigma = 2.2';
in.datasetbasename = 'Kin8nmcompactsigma2pt2';
in.datasetdir = datadirname;
in.methods = allmethods;

generate_dataset(in);
clear in;

%% Spam

load 'spam-distance-matrix'
in.sigma = 3;
in.d = 57;
in.cutoff = cutoffmultiplier*in.sigma;

in.A = generate_compact_RBF_kernel(spam_dist, in.sigma, in.d, in.cutoff);
in.sparsity = nnz(in.A)/prod(size(in.A));
clear spam_dist;
in.linearkernelflag = 0;
in.k = 10;
in.chunk = chunk;
in.vareps = vareps;
in.tol = tol;
in.lvals = genlvals(in.k, 150);
in.p = p;
in.q = q;

in.descr = 'Spam compact RBF kernel with sigma = 3';
in.datasetbasename = 'Spamcompactsigma3';
in.datasetdir = datadirname;
in.methods = allmethods;

generate_dataset(in);
clear in;

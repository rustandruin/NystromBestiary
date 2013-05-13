%% Runs randEig sketches on several datasets and stores 
% results

% specify which methods to use
allmethods = {'simple', 'srft', 'gaussian', 'levscore', 'froblev', ...
    'speclev', 'approxlev', ...
    'eig_simple', 'eig_srft', 'eig_gaussian', 'eig_levscore', ...
    'eig_froblev', 'eig_speclev', 'eig_approxlev', 'eignys_simple', ...
    'eignys_srft', 'eignys_gaussian', 'eignys_levscore', ...
    'eignys_froblev', 'eignys_speclev', 'eignys_approxlev'};
% specify the directory to output the results into
datadirname = 'densified';

q = 30; %6 % number of experiments to run
numpts = 12; %12 % number of different column samples to take
chunk = 10; % how often to reorthogonalize power method
vareps = 1/3; % ignored, is the epsilon parameter in some methods
tol = .01; % the tolerance to seek the leverage scores to in the power method
cutoffmultiplier = 3; % parameter for the generation of compact rbf matrices

% generates the set of numbers of column samples to be taken
genlvals = @(startl, endl) ceil(linspace(startl, endl, numpts));

%% Linear
% Dexter

load 'dextertestdata.sparse'
in.A = normalize_kernel_data(spconvert(dextertestdata));
clear dextertestdata;
in.linearkernelflag = 1;
in.k = 8;
in.chunk = chunk;
in.vareps = vareps;
in.tol = tol;
in.lvals = genlvals(in.k, 100);
in.q = q;

in.descr = '(whitened) Dexter linear kernel';
in.datasetbasename = 'Dexterrank8-eig';
in.datasetdir = datadirname;
in.methods = allmethods;

generate_dataset(in);
clear in;

%% RBF
% Abalone, sigma = .15

in.sigma = .15;
load 'abalone_distance_matrix'; 
in.A = generate_RBF_kernel(X, in.sigma);
clear X;

in.linearkernelflag = 0;
in.k = 20;
in.chunk = chunk;
in.vareps = vareps;
in.tol = tol;
in.lvals = genlvals(in.k, 150);
in.q = q;

in.descr = 'Abalone RBF kernel with sigma = .15';
in.datasetbasename = 'Abalonesigmapt15-eig';
in.datasetdir = datadirname;
in.methods = allmethods;

generate_dataset(in);
clear in;

%% Compact RBF
% Wine, sigma = 1

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
in.q = q;

in.descr = 'Wine compact RBF kernel with sigma = 1';
in.datasetbasename = 'Winecompactsigma1-eig';
in.datasetdir = datadirname;
in.methods = allmethods;

generate_dataset(in);
clear in;

%% Laplacian
% Gnutella rank 20

[~,~,in.A] = read_snap_data('p2p-Gnutella06.txt');
in.linearkernelflag = 0;
in.k = 20;
in.chunk = chunk;
in.vareps = vareps;
in.tol = tol;
in.lvals = genlvals(in.k, 180);
in.q = q;

in.descr = 'Gnutella normalized graph Laplacian';
in.datasetbasename = 'Gnutellarank20-eig';
in.datasetdir = datadirname;
in.methods = allmethods;

generate_dataset(in);
clear in;


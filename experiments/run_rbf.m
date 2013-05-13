%% Runs Nystrom extensions on various graph laplacian datasets and stores 
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

% generates the set of numbers of column samples to be taken
genlvals = @(startl, endl) ceil(linspace(startl, endl, numpts));


%% Abalone, sigma = 1

in.sigma = 1;
load 'abalone_distance_matrix'; 
in.A = generate_RBF_kernel(X, in.sigma);
clear X;

in.linearkernelflag = 0;
in.k = 20;
in.chunk = chunk;
in.vareps = vareps;
in.tol = tol;
in.lvals = genlvals(in.k, 150);
in.p = p;
in.q = q;

in.descr = 'Abalone RBF kernel with sigma = 1';
in.datasetbasename = 'Abalonesigma1';
in.datasetdir = datadirname;
in.methods = allmethods;

generate_dataset(in);
clear in;

%% Abalone, sigma = .15
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
in.p = p;
in.q = q;

in.descr = 'Abalone RBF kernel with sigma = .15';
in.datasetbasename = 'Abalonesigmapt15';
in.datasetdir = datadirname;
in.methods = allmethods;

generate_dataset(in);
clear in;

%% Wine, sigma = 1

load 'winequality_distance_matrix';
in.sigma = 1;
in.A = generate_RBF_kernel(X, in.sigma);
clear X;

in.linearkernelflag = 0;
in.k = 20;
in.chunk = chunk;
in.vareps = vareps;
in.tol = tol;
in.lvals = genlvals(in.k, 150);
in.p = p;
in.q = q;

in.descr = 'Wine RBF kernel with sigma = 1';
in.datasetbasename = 'Winesigma1';
in.datasetdir = datadirname;
in.methods = allmethods;

generate_dataset(in);
clear in;

%% Wine, sigma = 2.1

load 'winequality_distance_matrix';
in.sigma = 2.1;
in.A = generate_RBF_kernel(X, in.sigma);
clear X;

in.linearkernelflag = 0;
in.k = 20;
in.chunk = chunk;
in.vareps = vareps;
in.tol = tol;
in.lvals = genlvals(in.k, 150);
in.p = p;
in.q = q;

in.descr = 'Wine RBF kernel with sigma = 2.1';
in.datasetbasename = 'Winesigma2pt1';
in.datasetdir = datadirname;
in.methods = allmethods;

generate_dataset(in);
clear in;

%% Kin8nm

load 'kin8nm_distance_matrix'
in.sigma = 2.2;
in.A = generate_RBF_kernel(kin8nm_dist, in.sigma);
clear kin8nm_dist;

in.linearkernelflag = 0;
in.k = 40;
in.chunk = chunk;
in.vareps = vareps;
in.tol = tol;
in.lvals = genlvals(in.k, 150);
in.p = p;
in.q = q;

in.descr = 'Kin8nm RBF kernel with sigma = 2.1';
in.datasetbasename = 'Kin8nmsigma2pt1';
in.datasetdir = datadirname;
in.methods = allmethods;

generate_dataset(in);
clear in;

%% Spam

load 'spam-distance-matrix'
dataset = 'spamsigma3';
in.sigma = 3;
in.A = generate_RBF_kernel(spam_dist, in.sigma);
clear spam_dist;

in.linearkernelflag = 0;
in.k = 10;
in.chunk = chunk;
in.vareps = vareps;
in.tol = tol;
in.lvals = genlvals(in.k, 150);
in.p = p;
in.q = q;

in.descr = 'Spam RBF kernel with sigma = 3';
in.datasetbasename = 'Spamsigma2pt1';
in.datasetdir = datadirname;
in.methods = allmethods;

generate_dataset(in);
clear in;

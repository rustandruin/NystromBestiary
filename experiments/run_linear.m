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

%% Dexter

load 'dextertestdata.sparse'
in.A = normalize_kernel_data(spconvert(dextertestdata));
clear dextertestdata;
in.linearkernelflag = 1;
in.k = 8;
in.chunk = chunk;
in.vareps = vareps;
in.tol = tol;
in.lvals = genlvals(in.k, 100);
in.p = p;
in.q = q;

in.descr = '(whitened) Dexter linear kernel';
in.datasetbasename = 'Dexterrank8';
in.datasetdir = datadirname;
in.methods = allmethods;

generate_dataset(in);
clear in;

%% Protein

load 'protein_data';
in.A = normalize_kernel_data(X);
clear X;
in.linearkernelflag = 1;
in.k = 10; 
in.chunk = chunk;
in.vareps = vareps;
in.tol = tol;
in.lvals = genlvals(in.k, 110);
in.p = p;
in.q = q;

in.descr = '(whitened) Protein linear kernel';
in.datasetbasename = 'Proteinrank10';
in.datasetdir = datadirname;
in.methods = allmethods;

generate_dataset(in);
clear in;

%% SNPS 
load 'snpsdata'
A = A(1:(end-3), 1:size(A,2)-3); % The last three columns and rows are all 0s
in.A = normalize_kernel_data(A);
clear A;
in.linearkernelflag = 1;
in.k = 5; 
in.chunk = chunk;
in.vareps = vareps;
in.tol = tol;
in.lvals = genlvals(in.k, 46);
in.p = p;
in.q = q;

in.descr = '(whitened) SNPS linear kernel';
in.datasetbasename = 'SNPSrank5';
in.datasetdir = datadirname;
in.methods = allmethods;

generate_dataset(in);
clear in;

%% Gisette

load 'gisette_train.data';
A = normalize_kernel_data(gisette_train);
in.A = A*A';
clear gisette_train A;
in.linearkernelflag = 0;
in.k = 12; 
in.chunk = chunk;
in.vareps = vareps;
in.tol = tol;
in.lvals = genlvals(in.k, 150);
in.p = p;
in.q = q;

in.descr = '(whitened) Gisette linear kernel';
in.datasetbasename = 'Gisetterank12';
in.datasetdir = datadirname;
in.methods = allmethods;

generate_dataset(in);
clear in;

%% Cranfield Term document matrix

load 'cranfield-term-doc'
in.A = normalize_kernel_data(A_cran);
clear A_cran;
in.linearkernelflag = 1;
in.k = 3; 
in.chunk = chunk;
in.vareps = vareps;
in.tol = tol;
in.lvals = genlvals(in.k, 50);
in.p = p;
in.q = q;

in.descr = '(whitened) Cranfield linear kernel';
in.datasetbasename = 'Cranfieldrank3';
in.datasetdir = datadirname;
in.methods = allmethods;

generate_dataset(in);
clear in;

%% Medline Term Document Matrix

load 'medline-term-doc'
in.A = A_med;
clear A_med;
in.linearkernelflag = 1;
in.k = 15; 
in.chunk = chunk;
in.vareps = vareps;
in.tol = tol;
in.lvals = genlvals(in.k, 60);
in.p = p;
in.q = q;

in.descr = '(nonwhitened) Medline linear kernel';
in.datasetbasename = 'Medlinerank15';
in.datasetdir = datadirname;
in.methods = allmethods;

generate_dataset(in);
clear in;

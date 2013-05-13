%% Runs Nystrom extensions on various graph laplacian datasets and stores 
% results

% specify which methods to use
allmethods = {'simple', 'srft', 'gaussian', 'levscore', 'froblev', ...
    'speclev', 'approxlev'};
% specify the directory to output the results into
datadirname = 'outputs';

p = .9; % for the mixedprobs method
q = 6; %6 % number of experiments to run
numpts = 12; %12 % number of different column samples to take
chunk = 10; % how often to reorthogonalize power method
vareps = 1/3; % ignored, is the epsilon parameter in some methods
tol = .01; % the tolerance to seek the leverage scores to in the power method
cutoffmultiplier = 3; % parameter for the generation of compact rbf matrices

% generates the set of numbers of column samples to be taken
genlvals = @(startl, endl) ceil(linspace(startl, endl, numpts));

%% GR rank 60
[~,~,in.A] = read_snap_data('CA-GrQc.txt');
in.linearkernelflag = 0;
in.k = 60;
in.chunk = chunk;
in.vareps = vareps;
in.tol = tol;
in.lvals = genlvals(in.k, 180);
in.p = p;
in.q = q;

in.descr = 'GR normalized graph Laplacian';
in.datasetbasename = 'GRrank60';
in.datasetdir = datadirname;
in.methods = allmethods;

generate_dataset(in);
clear in;

%% GR rank 20
[~,~,in.A] = read_snap_data('CA-GrQc.txt');
in.linearkernelflag = 0;
in.k = 20;
in.chunk = chunk;
in.vareps = vareps;
in.tol = tol;
in.lvals = genlvals(in.k, 180);
in.p = p;
in.q = q;

in.descr = 'GR normalized graph Laplacian';
in.datasetbasename = 'GRrank20';
in.datasetdir = datadirname;
in.methods = allmethods;

generate_dataset(in);
clear in;

%% HEP rank 60
[~,~,in.A] = read_snap_data('CA-HepTh.txt');
in.linearkernelflag = 0;
in.k = 60;
in.chunk = chunk;
in.vareps = vareps;
in.tol = tol;
in.lvals = genlvals(in.k, 180);
in.p = p;
in.q = q;

in.descr = 'HEP normalized graph Laplacian';
in.datasetbasename = 'HEPrank60';
in.datasetdir = datadirname;
in.methods = allmethods;

generate_dataset(in);
clear in;

% HEP rank 20   
[~,~,in.A] = read_snap_data('CA-HepTh.txt');
in.linearkernelflag = 0;
in.k = 20;
in.chunk = chunk;
in.vareps = vareps;
in.tol = tol;
in.lvals = genlvals(in.k, 180);
in.p = p;
in.q = q;

in.descr = 'HEP normalized graph Laplacian';
in.datasetbasename = 'HEPrank20';
in.datasetdir = datadirname;
in.methods = allmethods;

generate_dataset(in);
clear in;

%% Enron 

% set up the matrix A
load 'email-Enron.mat'
sparseA = Problem.A;
n = size(sparseA, 1);
clear Problem;
degrees = sum(sparseA);
Dnegsqrt = diag(degrees.^(-1/2));
A = sparse(1:n, 1:n, 1) - Dnegsqrt*sparseA*Dnegsqrt; % the normalized Laplacian
clear n degrees Dnegsqrt
% make A smaller so the SRFT can be taken in memory
in.A = A(1:10e3, 1:10e3);
clear A sparseA;

% Enron rank 60

in.linearkernelflag = 0;
in.k = 60;
in.chunk = chunk;
in.vareps = vareps;
in.tol = tol;
in.lvals = genlvals(in.k, 180);
in.p = p;
in.q = q;

in.descr = 'Enron normalized graph Laplacian';
in.datasetbasename = 'Enronrank60';
in.datasetdir = datadirname;
in.methods = allmethods;

generate_dataset(in);
A = in.A;
clear in;

% Enron rank 20
in.A = A;
in.linearkernelflag = 0;
in.k = 20;
in.chunk = chunk;
in.vareps = vareps;
in.tol = tol;
in.lvals = genlvals(in.k, 180);
in.p = p;
in.q = q;

in.descr = 'Enron normalized graph Laplacian';
in.datasetbasename = 'Enronrank20';
in.datasetdir = datadirname;
in.methods = allmethods;

generate_dataset(in);
clear in;

%% Gnutella rank 60

[~,~,in.A] = read_snap_data('p2p-Gnutella06.txt');
in.linearkernelflag = 0;
in.k = 60;
in.chunk = chunk;
in.vareps = vareps;
in.tol = tol;
in.lvals = genlvals(in.k, 180);
in.p = p;
in.q = q;

in.descr = 'Gnutella normalized graph Laplacian';
in.datasetbasename = 'Gnutellarank60';
in.datasetdir = datadirname;
in.methods = allmethods;

generate_dataset(in);
clear in;

%% Gnutella rank 20
[~,~,in.A] = read_snap_data('p2p-Gnutella06.txt');
in.linearkernelflag = 0;
in.k = 20;
in.chunk = chunk;
in.vareps = vareps;
in.tol = tol;
in.lvals = genlvals(in.k, 180);
in.p = p;
in.q = q;

in.descr = 'Gnutella normalized graph Laplacian';
in.datasetbasename = 'Gnutellarank20';
in.datasetdir = datadirname;
in.methods = allmethods;

generate_dataset(in);
clear in;

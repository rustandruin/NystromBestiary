%% l = 6*k*log k for all datasets
genlvals = @(k) ceil(6*k*log(k));
q = 10;
datadirname = 'theory_practice_comparison_results_feb27';
comparisonmethods = {'simple', 'srft', 'gaussian', 'levscore'};

%% Enron, k=60

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

in.linearkernelflag = 0;
in.k = 60;
in.lvals = genlvals(in.k);
in.q = q;

in.descr = 'Enron normalized graph Laplacian';
in.datasetbasename = 'Enronrank60';
in.datasetdir = datadirname;
in.methods = comparisonmethods;

enrondata = generate_dataset(in);
clear in;

%% Protein, k = 10

load 'protein_data';
in.A = normalize_kernel_data(X);
clear X;
in.linearkernelflag = 1;
in.k = 10; % captures 42% of the frob norm, doubling to 20 captures only 51%
in.lvals = genlvals(in.k);
in.q = q;

in.descr = '(whitened) Protein linear kernel';
in.datasetbasename = 'Proteinrank10';
in.datasetdir = datadirname;
in.methods = comparisonmethods;

proteindata = generate_dataset(in);
clear in;

%% AbaloneD, sigma = .15, k = 20

in.sigma = .15;
load 'abalone_distance_matrix'; 
in.A = generate_RBF_kernel(X, in.sigma);
clear X;
in.linearkernelflag = 0;
in.k = 20;
in.lvals = genlvals(in.k);
in.q = q;

in.descr = 'Abalone RBF kernel with sigma = 1';
in.datasetbasename = 'Abalonesigma1';
in.datasetdir = datadirname;
in.methods = comparisonmethods;

abaloneddata = generate_dataset(in);
clear in;

%% WineS, sigma = 1, k = 20
cutoffmultiplier = 3;

load 'winequality_distance_matrix';
in.sigma = 1; 
in.d = 12;
in.cutoff = cutoffmultiplier*in.sigma;
in.A = generate_compact_RBF_kernel(X, in.sigma, in.d, in.cutoff);
in.sparsity = nnz(in.A)/prod(size(in.A));
clear X;
in.linearkernelflag = 0;
in.k = 20;
in.lvals = genlvals(in.k);
in.q = q;

in.descr = 'Wine compact RBF kernel with sigma = 1';
in.datasetbasename = 'Winecompactsigma1';
in.datasetdir = datadirname;
in.methods = comparisonmethods;

winesdata = generate_dataset(in);
clear in;

%% Save the data generated up to this point
save(fullfile(datadirname, 'entiredataset.mat'), ...
     'enrondata', 'proteindata', 'abaloneddata', 'winesdata');
 
%% Output the predicted errors and the actual errors

% Fix delta
delta = 1/2;

% Enron

% parameters for the Talwalkar bounds for this dataset
enronKmax = 1; % largest entry on diagonal of kernel
enrondmax = 2; % largest sqrt(A_{ii} + A_{jj} - 2A_{ij}) over i,j

bounds = predict_errors(enrondata, delta, enronKmax, enrondmax);
fprintf('\nDataset: Enron, l=%d\n', enrondata.in.lvals);
fprintf('Drineas: B2 = %f,\tBF = %f\n', ...
    bounds.drineas2bound, bounds.drineasFbound);
fprintf('Wolfe: BT = %f\n', ...
    bounds.belabbastrbound);
fprintf('Talwalkar: B2 = %f,\tBF = %f\n', ...
    bounds.talwalkar2bound, bounds.talwalkarFbound);
fprintf('Lemma 1: eps = %.2f,\tB2 = %f,\tBF = %f,\tBT = %f\n', ...
    bounds.levepsilon, bounds.lev2bound, ...
    bounds.levFbound, bounds.levtrbound);
fprintf('Lemma 2: eps = %.2f,\tB2 = %f,\tBF = %f,\tBT = %f\n', ...
    bounds.srftepsilon, bounds.srft2bound, ...
    bounds.srftFbound, bounds.srfttrbound);
fprintf('Lemma 3: eps = %.2f,\tB2 = %f,\tBF = %f,\tBT = %f\n', ...
    bounds.gaussepsilon, bounds.gauss2bound, ...
    bounds.gaussFbound, bounds.gausstrbound);
fprintf('Lemma 4: eps = %.2f,\tB2 = %f,\tBF = %f,\tBT = %f\n', ...
    bounds.nystromepsilon, bounds.nystrom2bound, ...
    bounds.nystromFbound, bounds.nystromtrbound);

% Protein

proteinKmax = 1;
proteindmax = 1.538;

bounds = predict_errors(proteindata, delta, proteinKmax, proteindmax);
fprintf('\nDataset: Protein\n');
fprintf('Drineas: B2 = %f,\tBF = %f\n', ...
    bounds.drineas2bound, bounds.drineasFbound);
fprintf('Wolfe: BT = %f\n', ...
    bounds.belabbastrbound);
fprintf('Talwalkar: B2 = %f,\tBF = %f\n', ...
    bounds.talwalkar2bound, bounds.talwalkarFbound);
fprintf('Lemma 1: eps = %.2f,\tB2 = %f,\tBF = %f,\tBT = %f\n', ...
    bounds.levepsilon, bounds.lev2bound, ...
    bounds.levFbound, bounds.levtrbound);
fprintf('Lemma 2: eps = %.2f,\tB2 = %f,\tBF = %f,\tBT = %f\n', ...
    bounds.srftepsilon, bounds.srft2bound, ...
    bounds.srftFbound, bounds.srfttrbound);
fprintf('Lemma 3: eps = %.2f,\tB2 = %f,\tBF = %f,\tBT = %f\n', ...
    bounds.gaussepsilon, bounds.gauss2bound, ...
    bounds.gaussFbound, bounds.gausstrbound);
fprintf('Lemma 4: eps = %.2f,\tB2 = %f,\tBF = %f,\tBT = %f\n', ...
    bounds.nystromepsilon, bounds.nystrom2bound, ...
    bounds.nystromFbound, bounds.nystromtrbound);

% AbaloneD

abaloneDKmax = 1;
abaloneDdmax = 1.538;

bounds = predict_errors(abaloneddata, delta, abaloneDKmax, abaloneDdmax);
fprintf('\nDataset: AbaloneD\n');
fprintf('Drineas: B2 = %f,\tBF = %f\n', ...
    bounds.drineas2bound, bounds.drineasFbound);
fprintf('Wolfe: BT = %f\n', ...
    bounds.belabbastrbound);
fprintf('Talwalkar: B2 = %f,\tBF = %f\n', ...
    bounds.talwalkar2bound, bounds.talwalkarFbound);
fprintf('Lemma 1: eps = %.2f,\tB2 = %f,\tBF = %f,\tBT = %f\n', ...
    bounds.levepsilon, bounds.lev2bound, ...
    bounds.levFbound, bounds.levtrbound);
fprintf('Lemma 2: eps = %.2f,\tB2 = %f,\tBF = %f,\tBT = %f\n', ...
    bounds.srftepsilon, bounds.srft2bound, ...
    bounds.srftFbound, bounds.srfttrbound);
fprintf('Lemma 3: eps = %.2f,\tB2 = %f,\tBF = %f,\tBT = %f\n', ...
    bounds.gaussepsilon, bounds.gauss2bound, ...
    bounds.gaussFbound, bounds.gausstrbound);
fprintf('Lemma 4: eps = %.2f,\tB2 = %f,\tBF = %f,\tBT = %f\n', ...
    bounds.nystromepsilon, bounds.nystrom2bound, ...
    bounds.nystromFbound, bounds.nystromtrbound);

% WineS

wineSKmax = 1;
wineSdmax = sqrt(2);

bounds = predict_errors(winesdata, delta, wineSKmax, wineSdmax);
fprintf('\nDataset: WineS\n');
fprintf('Drineas: B2 = %f,\tBF = %f\n', ...
    bounds.drineas2bound, bounds.drineasFbound);
fprintf('Wolfe: BT = %f\n', ...
    bounds.belabbastrbound);
fprintf('Talwalkar: B2 = %f,\tBF = %f\n', ...
    bounds.talwalkar2bound, bounds.talwalkarFbound);
fprintf('Lemma 1: eps = %.2f,\tB2 = %f,\tBF = %f,\tBT = %f\n', ...
    bounds.levepsilon, bounds.lev2bound, ...
    bounds.levFbound, bounds.levtrbound);
fprintf('Lemma 2: eps = %.2f,\tB2 = %f,\tBF = %f,\tBT = %f\n', ...
    bounds.srftepsilon, bounds.srft2bound, ...
    bounds.srftFbound, bounds.srfttrbound);
fprintf('Lemma 3: eps = %.2f,\tB2 = %f,\tBF = %f,\tBT = %f\n', ...
    bounds.gaussepsilon, bounds.gauss2bound, ...
    bounds.gaussFbound, bounds.gausstrbound);
fprintf('Lemma 4: eps = %.2f,\tB2 = %f,\tBF = %f,\tBT = %f\n', ...
    bounds.nystromepsilon, bounds.nystrom2bound, ...
    bounds.nystromFbound, bounds.nystromtrbound);

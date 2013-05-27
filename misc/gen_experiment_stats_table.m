%% l = 6*k*log k for all datasets
genlvals = @(k,n) [k + 8, ceil(k*log(k)), ceil(k*log(n))];
q = 30;
datadirname = 'theory_practice_comparison_results_may22';
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
clear degrees Dnegsqrt
% make A smaller so the SRFT can be taken in memory
in.A = A(1:10e3, 1:10e3);
clear A sparseA;

in.linearkernelflag = 0;
in.k = 60;
in.lvals = genlvals(in.k,n);
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
n = size(in.A,1);
clear X;
in.linearkernelflag = 1;
in.k = 10; % captures 42% of the frob norm, doubling to 20 captures only 51%
in.lvals = genlvals(in.k,n);
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
n = size(in.A, 1);
clear X;
in.linearkernelflag = 0;
in.k = 20;
in.lvals = genlvals(in.k,n);
in.q = q;

in.descr = 'Abalone RBF kernel with sigma = .15';
in.datasetbasename = 'Abalonesigmapt15';
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
n = size(in.A, 2);
in.sparsity = nnz(in.A)/prod(size(in.A));
clear X;
in.linearkernelflag = 0;
in.k = 20;
in.lvals = genlvals(in.k,n);
in.q = q;

in.descr = 'Wine compact RBF kernel with sigma = 1';
in.datasetbasename = 'Winecompactsigma1';
in.datasetdir = datadirname;
in.methods = comparisonmethods;

winesdata = generate_dataset(in);
clear in;

%% Save the data generated up to this point
save(fullfile(datadirname, 'entire-stats-dataset.mat'), ...
     'enrondata', 'proteindata', 'abaloneddata', 'winesdata');
 
%% Output the statistics of the errors


for dataset = {enrondata, proteindata, abaloneddata, winesdata}
    dataset = dataset{1};
    
    fprintf('\n\n%s data set\n', dataset.in.datasetbasename);
    fprintf('----------------------------------------\n\n');
    [nystrom, srft, gaussian, levscore] = genstats(dataset);
    
    fprintf('Nystrom, nonrestricted\n');
    fprintf('----------------------\n');
    fprintf('Spectral (min, mean, max)\t%.2f/%.2f/%.2f\t%.2f/%.2f/%.2f\t%.2f/%.2f/%.2f\n', nystrom.nonrestricted_spec_stats');
    fprintf('Frobenius (min, mean, max)\t%.2f/%.2f/%.2f\t%.2f/%.2f/%.2f\t%.2f/%.2f/%.2f\n', nystrom.nonrestricted_fro_stats');
    fprintf('Trace (min, mean, max)\t\t%.2f/%.2f/%.2f\t%.2f/%.2f/%.2f\t%.2f/%.2f/%.2f\n\n', nystrom.nonrestricted_tr_stats');
    
    fprintf('Nystrom, restricted\n');
    fprintf('-------------------\n');
    fprintf('Spectral (min, mean, max)\t%.2f/%.2f/%.2f\t%.2f/%.2f/%.2f\t%.2f/%.2f/%.2f\n', nystrom.restricted_spec_stats');
    fprintf('Frobenius (min, mean, max)\t%.2f/%.2f/%.2f\t%.2f/%.2f/%.2f\t%.2f/%.2f/%.2f\n', nystrom.restricted_fro_stats');
    fprintf('Trace (min, mean, max)\t\t%.2f/%.2f/%.2f\t%.2f/%.2f/%.2f\t%.2f/%.2f/%.2f\n\n', nystrom.restricted_tr_stats');
    
    fprintf('SRFT, nonrestricted\n');
    fprintf('-------------------\n');
    fprintf('Spectral (min, mean, max)\t%.2f/%.2f/%.2f\t%.2f/%.2f/%.2f\t%.2f/%.2f/%.2f\n', srft.nonrestricted_spec_stats');
    fprintf('Frobenius (min, mean, max)\t%.2f/%.2f/%.2f\t%.2f/%.2f/%.2f\t%.2f/%.2f/%.2f\n', srft.nonrestricted_fro_stats');
    fprintf('Trace (min, mean, max)\t\t%.2f/%.2f/%.2f\t%.2f/%.2f/%.2f\t%.2f/%.2f/%.2f\n\n', srft.nonrestricted_tr_stats');
    
    fprintf('SRFT, restricted\n');
    fprintf('----------------\n');
    fprintf('Spectral (min, mean, max)\t%.2f/%.2f/%.2f\t%.2f/%.2f/%.2f\t%.2f/%.2f/%.2f\n', srft.restricted_spec_stats');
    fprintf('Frobenius (min, mean, max)\t%.2f/%.2f/%.2f\t%.2f/%.2f/%.2f\t%.2f/%.2f/%.2f\n', srft.restricted_fro_stats');
    fprintf('Trace (min, mean, max)\t\t%.2f/%.2f/%.2f\t%.2f/%.2f/%.2f\t%.2f/%.2f/%.2f\n\n', srft.restricted_tr_stats');
    
    fprintf('Gaussian, nonrestricted\n');
    fprintf('-----------------------\n');
    fprintf('Spectral (min, mean, max)\t%.2f/%.2f/%.2f\t%.2f/%.2f/%.2f\t%.2f/%.2f/%.2f\n', gaussian.nonrestricted_spec_stats');
    fprintf('Frobenius (min, mean, max)\t%.2f/%.2f/%.2f\t%.2f/%.2f/%.2f\t%.2f/%.2f/%.2f\n', gaussian.nonrestricted_fro_stats');
    fprintf('Trace (min, mean, max)\t\t%.2f/%.2f/%.2f\t%.2f/%.2f/%.2f\t%.2f/%.2f/%.2f\n\n', gaussian.nonrestricted_tr_stats');
    
    fprintf('Gaussian, restricted\n');
    fprintf('--------------------\n');
    fprintf('Spectral (min, mean, max)\t%.2f/%.2f/%.2f\t%.2f/%.2f/%.2f\t%.2f/%.2f/%.2f\n', gaussian.restricted_spec_stats');
    fprintf('Frobenius (min, mean, max)\t%.2f/%.2f/%.2f\t%.2f/%.2f/%.2f\t%.2f/%.2f/%.2f\n', gaussian.restricted_fro_stats');
    fprintf('Trace (min, mean, max)\t\t%.2f/%.2f/%.2f\t%.2f/%.2f/%.2f\t%.2f/%.2f/%.2f\n\n', gaussian.restricted_tr_stats');
    
    fprintf('Levscore, nonrestricted\n');
    fprintf('-----------------------\n');
    fprintf('Spectral (min, mean, max)\t%.2f/%.2f/%.2f\t%.2f/%.2f/%.2f\t%.2f/%.2f/%.2f\n', levscore.nonrestricted_spec_stats');
    fprintf('Frobenius (min, mean, max)\t%.2f/%.2f/%.2f\t%.2f/%.2f/%.2f\t%.2f/%.2f/%.2f\n', levscore.nonrestricted_fro_stats');
    fprintf('Trace (min, mean, max)\t\t%.2f/%.2f/%.2f\t%.2f/%.2f/%.2f\t%.2f/%.2f/%.2f\n\n', levscore.nonrestricted_tr_stats');
    
    fprintf('Levscore, restricted\n');
    fprintf('--------------------\n');
    fprintf('Spectral (min, mean, max)\t%.2f/%.2f/%.2f\t%.2f/%.2f/%.2f\t%.2f/%.2f/%.2f\n', levscore.restricted_spec_stats');
    fprintf('Frobenius (min, mean, max)\t%.2f/%.2f/%.2f\t%.2f/%.2f/%.2f\t%.2f/%.2f/%.2f\n', levscore.restricted_fro_stats');
    fprintf('Trace (min, mean, max)\t\t%.2f/%.2f/%.2f\t%.2f/%.2f/%.2f\t%.2f/%.2f/%.2f\n', levscore.restricted_tr_stats');
    
end

%% generate LaTeX tables

for dataset = {enrondata, proteindata, abaloneddata, winesdata}
    dataset = dataset{1};
    
    fprintf('\n\n%s data set\n', dataset.in.datasetbasename);
    fprintf('----------------------------------------\n\n');
    [nystrom, srft, gaussian, levscore] = genstats(dataset);
    
    fprintf('\\hline\n');
    fprintf('$\\|\\matA - \\matC\\matW^\\pinv\\matC^\\transp\\|_2/\\|\\matA - \\matA_k\\|_2$ \\\\\n');
    fprintf('\\hline\n');
    fprintf('\\begin{tabular}{lccc}\n')
    fprintf('& $\\ell = k+8$ & $\\ell = k\\ln k$ & $\\ell = k \\ln n$ \\\\\n');
    fprintf('Nystr\\"om & %.3f/%.3f/%.3f & %.3f/%.3f/%.3f & %.3f/%.3f/%.3f \\\\\n', nystrom.nonrestricted_spec_stats');
    fprintf('SRFT sketch & %.3f/%.3f/%.3f & %.3f/%.3f/%.3f & %.3f/%.3f/%.3f \\\\\n', srft.nonrestricted_spec_stats');
    fprintf('Gaussian sketch & %.3f/%.3f/%.3f & %.3f/%.3f/%.3f & %.3f/%.3f/%.3f \\\\\n', gaussian.nonrestricted_spec_stats');
    fprintf('Leverage sketch & %.3f/%.3f/%.3f & %.3f/%.3f/%.3f & %.3f/%.3f/%.3f\n', levscore.nonrestricted_spec_stats');
    fprintf('\\end{tabular} \\\\\n');
    fprintf('\\hline\n');
    fprintf('$\\|\\matA - \\matC\\matW^\\pinv\\matC^\\transp\\|_{\\mathrm{F}}/\\|\\matA - \\matA_k\\|_{\\mathrm{F}}$ \\\\\n');
    fprintf('\\hline\n');
    fprintf('\\begin{tabular}{lccc}\n')
    fprintf('& $\\ell = k+8$ & $\\ell = k\\ln k$ & $\\ell = k \\ln n$ \\\\\n');
    fprintf('Nystr\\"om & %.3f/%.3f/%.3f & %.3f/%.3f/%.3f & %.3f/%.3f/%.3f \\\\\n', nystrom.nonrestricted_fro_stats');
    fprintf('SRFT sketch & %.3f/%.3f/%.3f & %.3f/%.3f/%.3f & %.3f/%.3f/%.3f \\\\\n', srft.nonrestricted_fro_stats');
    fprintf('Gaussian sketch & %.3f/%.3f/%.3f & %.3f/%.3f/%.3f & %.3f/%.3f/%.3f \\\\\n', gaussian.nonrestricted_fro_stats');
    fprintf('Leverage sketch & %.3f/%.3f/%.3f & %.3f/%.3f/%.3f & %.3f/%.3f/%.3f\n', levscore.nonrestricted_fro_stats');
    fprintf('\\end{tabular} \\\\ \n');
    fprintf('\\hline\n');
    fprintf('$\\|\\matA - \\matC\\matW^\\pinv\\matC^\\transp\\|_{\\star}/\\|\\matA - \\matA_k\\|_{\\star}$ \\\\\n');
    fprintf('\\hline\n');
    fprintf('\\begin{tabular}{lccc}\n')
    fprintf('& $\\ell = k+8$ & $\\ell = k\\ln k$ & $\\ell = k \\ln n$ \\\\\n');
    fprintf('Nystr\\"om & %.3f/%.3f/%.3f & %.3f/%.3f/%.3f & %.3f/%.3f/%.3f \\\\\n', nystrom.nonrestricted_tr_stats');
    fprintf('SRFT sketch & %.3f/%.3f/%.3f & %.3f/%.3f/%.3f & %.3f/%.3f/%.3f \\\\\n', srft.nonrestricted_tr_stats');
    fprintf('Gaussian sketch & %.3f/%.3f/%.3f & %.3f/%.3f/%.3f & %.3f/%.3f/%.3f \\\\\n', gaussian.nonrestricted_tr_stats');
    fprintf('Leverage sketch & %.3f/%.3f/%.3f & %.3f/%.3f/%.3f & %.3f/%.3f/%.3f\n', levscore.nonrestricted_tr_stats');
    fprintf('\\end{tabular} \\\\ \n');
    fprintf('\\hline\n');
end
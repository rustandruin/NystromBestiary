% Note: this will generate about 1000Mb of data

% Create the distance matrix for the Abalone dataset
clear
load abalone_dataset
X = generate_distance_matrix(abaloneInputs');
[basepath, ~,~] = fileparts(mfilename('fullpath'));
fname = fullfile(basepath, 'abalone_distance_matrix.mat');
save(fname);

% Create the distance matrix for the Wine dataset
clear
load winequality_data
X = generate_distance_matrix(X);
[basepath, ~,~] = fileparts(mfilename('fullpath'));
fname = fullfile(basepath, 'winequality_distance_matrix.mat');
save(fname);

% Create the distance matrix for the Kin8nm dataset
clear
load kin8nm.data
kin8nm_dist = generate_distance_matrix(kin8nm);
[basepath, ~,~] = fileparts(mfilename('fullpath'));
fname = fullfile(basepath, 'kin8nm_distance_matrix.mat');
save(fname);

% Create the distance matrix for the Spam dataset
clear 
load spam.data
spam_dist = generate_distance_matrix(spam);
[basepath, ~,~] = fileparts(mfilename('fullpath'));
fname = fullfile(basepath, 'spam-distance-matrix.mat');
save(fname);

clear




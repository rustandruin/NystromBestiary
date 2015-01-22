function savedata = generate_dataset(in)
%
% savedata = generate_dataset(in)
%
% Stores the error and timing information from running chosen SPSD
% sketching methods on a given dataset with given parameters using a 
% specified range of column samples
%
% in is a structure with (at least) the following fields:
% - descr, a string describing this dataset
% - datasetbasename, a base filename for this dataset
% - datasetdir, the directory to store this dataset in
% - methods, a cell array of the sketching methods to apply; see below for
%   valid methods
% - A, a matrix
% - linearkernelflag, is 0 if A is PSD, and 1 is AA^T is the PSD matrix to
%   be approximated
% - k, the target rank of the approximation
% - lvals, a vector specifying the numbers of column samples to use
% - q, the number of times to repeat each sketching method for each number of
%  column samples
%
% Other fields may be required to be present, depending on the sketching
% methods specified. Valid methods (case insensitive):
% 'simple' - uniform column sampling without replacement
% 'srft' - SRFT mixture-based
% 'gaussian' - Gaussian mixture-based
% 'levscore' - leverage score based column sampling
% 'froblev' - frobenius sketch approximate leverage score sampling
%  requires additional field vareps, a number between 0 and 1
% 'speclev' - spectral sketch approximate leverage score sampling
%  requires additional fields vareps and chunk, how many iterations to take
%  before reorthogonalization 
% 'approxlev' - power method approximate leverage score sampling
%  requires additional fields vareps and chunk
%  'tallthinlevscoreapprox' - uses Alg 1 from Mahoney to do approximate
%  leverage score-based column sampling
%
%  Any of the above methods prefixed with 'eig_' runs the randomized SVD
%  variant of the sketch, i.e. approximates A with QQ*AQQ*, where QR = AS
%  for some sketching matrix S. Any prefixed with 'eignys_' runs the 
%  randomized SVD-based Nystrom extension variant of the sketch, i.e.
%  approximates A with (AQ)pinv(Q*AQ)(AQ)*
%  
%
% The return struct savedata contains all the information saved from the
% skteching experiments run.
%

in.methods = lower(in.methods);
wantq = @(methodname) any(strcmp(methodname, in.methods));
savedata.in = in;
datasetfname = fullfile(in.datasetdir, in.datasetbasename);

% compute the leverage scores and optimal rank-k approximation errors
if in.linearkernelflag == 0 % A is a PSD matrix
    tic
       [U, Sigma] = orderedeigs(in.A, in.k+1);
       U1t = U(:, 1:in.k)';
       savedata.levscores = sum(U1t.*U1t);
    in.levscorecomputationtime = toc;
    in.levscoreprobs = savedata.levscores/in.k;
    savedata.topspectrum = diag(Sigma(1:in.k,1:in.k));
    savedata.optspecerr = Sigma(in.k+1,in.k+1);
    savedata.optfroerr = sqrt(norm(in.A, 'fro')^2 - sum(savedata.topspectrum.^2));
    savedata.opttrerr = trace(in.A) - sum(savedata.topspectrum);
else % the actual PSD matrix is AA^T
    tic
       [U, Sigma, ~] = svds(in.A, in.k+1);
       U1t = U(:, 1:in.k)';
       savedata.levscores = sum(U1t.*U1t);
    in.levscoreprobs = savedata.levscores/in.k;
    in.levscorecomputationtime = toc;
    savedata.topspectrum = diag(Sigma(1:in.k, 1:in.k)).^2;
    savedata.optspecerr = Sigma(in.k+1,in.k+1)^2;
    A = in.A*in.A';
    savedata.optfroerr = sqrt(norm(A, 'fro')^2 - sum(savedata.topspectrum.^2));
    savedata.opttrerr = trace(A) - sum(savedata.topspectrum);
end

% if we use the tall thin leverage score algorithm to compute Nystrom
% approximations, we should compare it to using the QR algorithm to compute
% the standard leverage score extension for fairness
if wantq('tallthinlevscoreapprox') || wantq('eig_tallthinlevscoreapprox')
    tic
        [Q,~] = qr(in.A,0);
        in.qrlevscoreprobs = sum(Q.*Q,2)/size(in.A,2);
    in.qrlevscorecomputationtime = toc;
end

% store the errors of the specified Nystrom methods
for lidx = 1:length(in.lvals)
    in.l = in.lvals(lidx);
    
    fprintf('Evaluating errors for %s, l = %d (%d of %d values)\n', in.datasetbasename, in.l, lidx, length(in.lvals));
    
    if wantq('testmethod')
        fprintf('...test method');
        testmethodData(lidx) = test_Nystrom(in);
    end
    
    if wantq('simple')
        fprintf('...simple\n');
        simpleData(lidx) = simple_Nystrom(in);
    end
    
    if wantq('srft')
        fprintf('...srft\n');
        srftData(lidx) = srft_Nystrom(in);
    end
    
    if wantq('gaussian')
        fprintf('...gaussian\n');
        gaussianData(lidx) = gaussian_Nystrom(in);
    end
    
    if wantq('levscore')
        fprintf('...levscore\n');
        levscoreData(lidx) = levscore_Nystrom(in);
    end
    
    if wantq('froblev')
        fprintf('...froblev\n');
        froblevData(lidx) = froblev_Nystrom(in);
    end
    
    if wantq('speclev')
        fprintf('...speclev\n');
        speclevData(lidx) = speclev_Nystrom(in);
    end
    
    if wantq('approxlev')
        fprintf('...approxlev\n');
        approxlevData(lidx) = approxlev_Nystrom(in);
    end
    
    if wantq('tallthinlevscoreapprox')
        fprintf('...tallthinlevscoreapprox\n');
            tallthinData(lidx) = tallthin_Nystrom(in);
        fprintf('...qrlevscore for comparison to tallthinalg\n');
            qrlevscoreData(lidx) = qrlevscore_Nystrom(in);
    end
    
    if wantq('testmethod')
        fprintf('...test method');
        testmethodData(lidx) = test_Nystrom(in);
    end
    
    if wantq('eig_simple')
        fprintf('...eig_simple\n');
        eig_simpleData(lidx) = simple_randEig(in);
    end
    
    if wantq('eig_srft')
        fprintf('...eig_srft\n');
        eig_srftData(lidx) = srft_randEig(in);
    end
    
    if wantq('eig_gaussian')
        fprintf('...eig_gaussian\n');
        eig_gaussianData(lidx) = gaussian_randEig(in);
    end
    
    if wantq('eig_levscore')
        fprintf('...eig_levscore\n');
        eig_levscoreData(lidx) = levscore_randEig(in);
    end
    
    if wantq('eig_froblev')
        fprintf('...eig_froblev\n');
        eig_froblevData(lidx) = froblev_randEig(in);
    end
    
    if wantq('eig_speclev')
        fprintf('...eig_speclev\n');
        eig_speclevData(lidx) = speclev_randEig(in);
    end
    
    if wantq('eig_approxlev')
        fprintf('...eig_approxlev\n');
        eig_approxlevData(lidx) = approxlev_randEig(in);
    end

    if wantq('eig_tallthinlevscoreapprox')
        fprintf('...eig_tallthinlevscoreapprox\n');
            eig_tallthinData(lidx) = tallthin_randEig(in);
        fprintf('...qrlevscore for comparison to tallthinalg\n');
            eig_qrlevscoreData(lidx) = qrlevscore_randEig(in);
    end
    
    if wantq('eignys_simple')
        fprintf('...eignys_simple\n');
        eignys_simpleData(lidx) = simple_randEig2(in);
    end
    
    if wantq('eignys_srft')
        fprintf('...eignys_srft\n');
        eignys_srftData(lidx) = srft_randEig2(in);
    end
    
    if wantq('eignys_gaussian')
        fprintf('...eignys_gaussian\n');
        eignys_gaussianData(lidx) = gaussian_randEig2(in);
    end
    
    if wantq('eignys_levscore')
        fprintf('...eignys_levscore\n');
        eignys_levscoreData(lidx) = levscore_randEig2(in);
    end
    
    if wantq('eignys_froblev')
        fprintf('...eignys_froblev\n');
        eignys_froblevData(lidx) = froblev_randEig2(in);
    end
    
    if wantq('eignys_speclev')
        fprintf('...eignys_speclev\n');
        eignys_speclevData(lidx) = speclev_randEig2(in);
    end
    
    if wantq('eignys_approxlev')
        fprintf('...eignys_approxlev\n');
        eignys_approxlevData(lidx) = approxlev_randEig2(in);
    end

    if wantq('eignys_tallthinlevscoreapprox')
        fprintf('...eignys_tallthinlevscoreapprox\n');
            eignys_tallthinData(lidx) = tallthin_randEig2(in);
        fprintf('...qrlevscore for comparison to tallthinalg\n');
            eignys_qrlevscoreData(lidx) = qrlevscore_randEig2(in);
    end
end

if wantq('testmethod')
    savedata.testmethodData = testmethodData;
end

if wantq('simple')
    savedata.simpleData = simpleData;
end

if wantq('srft')
    savedata.srftData = srftData;
end

if wantq('gaussian')
    savedata.gaussianData = gaussianData;
end

if wantq('levscore')
    savedata.levscoreData = levscoreData;
end

if wantq('froblev')
    savedata.froblevData = froblevData;
end

if wantq('speclev')
    savedata.speclevData = speclevData;
end

if wantq('approxlev')
    savedata.approxlevData = approxlevData;
end

if wantq('tallthinlevscoreapprox')
    savedata.tallthinData = tallthinData;
    savedata.qrlevscoreData = qrlevscoreData;
end

if wantq('eig_simple')
    savedata.eig_simpleData = eig_simpleData;
end

if wantq('eig_srft')
    savedata.eig_srftData = eig_srftData;
end

if wantq('eig_gaussian')
    savedata.eig_gaussianData = eig_gaussianData;
end

if wantq('eig_levscore')
    savedata.eig_levscoreData = eig_levscoreData;
end

if wantq('eig_froblev')
    savedata.eig_froblevData = eig_froblevData;
end

if wantq('eig_speclev')
    savedata.eig_speclevData = eig_speclevData;
end

if wantq('eig_approxlev')
    savedata.eig_approxlevData = eig_approxlevData;
end

if wantq('eig_tallthinlevscoreapprox')
    savedata.eig_tallthinData = eig_tallthinData;
    savedata.eig_qrlevscoreData = eig_qrlevscoreData;
end

if wantq('eignys_simple')
    savedata.eignys_simpleData = eignys_simpleData;
end

if wantq('eignys_srft')
    savedata.eignys_srftData = eignys_srftData;
end

if wantq('eignys_gaussian')
    savedata.eignys_gaussianData = eignys_gaussianData;
end

if wantq('eignys_levscore')
    savedata.eignys_levscoreData = eignys_levscoreData;
end

if wantq('eignys_froblev')
    savedata.eignys_froblevData = eignys_froblevData;
end

if wantq('eignys_speclev')
    savedata.eignys_speclevData = eignys_speclevData;
end

if wantq('eignys_approxlev')
    savedata.eignys_approxlevData = eignys_approxlevData;
end

if wantq('eignys_tallthinlevscoreapprox')
    savedata.eignys_tallthinData = eignys_tallthinData;
    savedata.eignys_qrlevscoreData = eignys_qrlevscoreData;
end

save(datasetfname, 'savedata');

end

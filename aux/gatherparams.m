% stats = gatherparams(fname)
%
% loads in a file of stored experiment data and returns a 
% structure with relevant statistics

function stats = gatherparams(fname)

load(fname);

k = savedata.in.k;

if ~savedata.in.linearkernelflag
    A = savedata.in.A;
else
    A = savedata.in.A * savedata.in.A';   
end

n = size(A,1);

topspectrum = svds(A, k+1);
frobnorm = norm(A, 'fro');
tracenorm = trace(A);
sortedlevscores = sort(savedata.levscores, 2, 'descend');

numericalrank = frobnorm^2/topspectrum(1)^2;
spectralratio = topspectrum(end)/topspectrum(end-1);
frobmass = 100*norm(topspectrum(1:k))/frobnorm;
tracemass = 100*sum(topspectrum(1:k))/tracenorm;
toplevscale = n/k*sortedlevscores(k);
normalizedstd = std(sortedlevscores)/mean(sortedlevscores);

thres = 1e-8; % consider leverage scores smaller than this to be 0
levscoredist = sortedlevscores(sortedlevscores >= thres)/k;
entdist = sum(levscoredist .* log2(1./levscoredist));

nnzpercentage = nnz(A)/prod(size(A))*100;

stats = struct('k', k, ...
               'nnzpercentage', nnzpercentage, ...
               'numericalrank', numericalrank, ...
               'spectralratio', spectralratio, ...
               'frobmass', frobmass, ...
               'tracemass', tracemass, ...
               'toplevscale', toplevscale, ...
               'normalizedstd', normalizedstd, ...
               'entdist', entdist);
end
    
    

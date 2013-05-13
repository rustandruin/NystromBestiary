function out = tallthin_Nystrom(in)
% out = tallthin_Nystrom(in)
%
% in is a structure with (at least) the following fields:
%  -A, an n-by-n PSD matrix, and
% -linearkernelflag, 1 if A=XX^T and X was passed in A's position, 0
%  otherwise, and
%  -l < n, the number of columns to sample and
%  -q >= 1 the number of times to repeat the experiment and
%  -chunk, how often to check for convergence of the approx leverage scores and
%  -tallthineps, the additive accuracy we desire in the leverage scores
%  -tallthinmethod, 'qr' or 'svd', indicating which method to use 
%
% out is a structure with the following fields:
%  -specerr, froerr, trerr: each matrices with two rows representing 
%  the spectral, frobenius, 
%  and trace norms of the errors in using q realization of an approximate
%  leverage score based Nystrom extension to A, using l columns sampled 
%  according to approximate leverage score probabilities. The first row corresponds to 
%  Nystrom extensions where the rank was not fixed, the second to Nystrom
%  extensions where the rank was fixed
% 
%  -timings, a two row vector of the time it took to run each experiment 
%  (i.e. form C,Winv,Wkinv), including the time to approximate the leverage 
%  scores
%
%  -numiters, a vector of the number of iterations it took to obtain the 
%   approximate leverage scores to the desired accuracy
%
%  -approxlevscores, a vector of one set of leverage score approximations
%
% See also APPROX_TALLTHIN_LEVSCORES

% % if testing
%out = dummy_extension(in);
%return;

n = size(in.A,1);
out.specerr = zeros(2,in.q);
out.froerr = zeros(2,in.q);
out.trerr = zeros(2,in.q);
out.timings = zeros(2,in.q);
out.numiters = zeros(1,in.q);

in.maxiters = ceil(log(in.vareps)/log(.9)); % This corresponds to guessing the spectral gap ratio is smaller than .9
out.maxiters = in.maxiters; % save this parameter

Id = speye(n);

if in.linearkernelflag
    At = in.A';
end

for iter=1:in.q
    tic
        % compute the approximate leverage scores
        [~,~,~,out.tallthinlevscores] = ...
            approx_tallthin_levscores(in.A, in.tallthineps, in.tallthinmethod);
        levscoreprobs = out.tallthinlevscores/size(in.A, 2); % the leverage scores are for the full rank (d)!

        % sample according to those leverage scores
        colindices = ones(1,in.l);
        for i=1:in.l
            colindices(i) = find(cumsum(levscoreprobs) >= rand(),1);
        end
        scalingfactors = levscoreprobs(colindices).^(1/2);
        S = Id(:,colindices)*diag(1./scalingfactors);

        if in.linearkernelflag == 0
            C = in.A*S;
        else
            C = in.A*(At*S);
        end

        W = full(S'*C);
    out.timings(:, iter) = toc;
    
    % time the eigenvalue decomposition of W
    tic
        [V,D] = orderedeig(W);
    decomptime = toc;
    
    % time to form Winv
    tic
        Winv = V*pinv(D)*V';
    out.timings(1, iter) = out.timings(1, iter) + decomptime + toc;
    
    % time to form Wkinv
    tic
        Wkinv = V(:, 1:in.k)*pinv(D(1:in.k, 1:in.k))*V(:,1:in.k)';
    out.timings(2, iter) = out.timings(2, iter) + decomptime + toc;
    
    % estimate the Nystrom approx erros
    [out.specerr(:,iter), out.froerr(:,iter), out.trerr(:,iter)] = ...
        estnorms(in, C, Winv, Wkinv);
end

end

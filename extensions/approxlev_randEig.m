function out = approxlev_randEig(in)
% out = approxlev_randEig(in)
%
% in is a structure with (at least) the following fields:
%  -A, an n-by-n PSD matrix, and
% -linearkernelflag, 1 if A=XX^T and X was passed in A's position, 0
%  otherwise, and
%  -k, the target rank, and 
%  -l < n, the number of columns to sample and
%  -q >= 1 the number of times to repeat the experiment and
%  -chunk, how often to check for convergence of the approx leverage scores and
%  -vareps, the additive accuracy we desire in the leverage scores
%
% out is a structure with the following fields:
%  -specerr, froerr, trerr: each matrices with two rows representing 
%  the spectral, frobenius, 
%  and trace norms of the errors in using q realization of the leverage
%  score based randEig sketch to A, using l columns sampled according to
%  approximate leverage score probabilities. The first row corresponds to 
%  randEig sketches where the rank was not fixed, the second to randEig sketches
%  where the rank was fixed
% 
%  -timings, a two row vector of the time it took to run each experiment 
%  (i.e. form C,Q,W,Wk), including the time to approximate the leverage 
%  scores
%
%  -numiters, a vector of the number of iterations it took to obtain the 
%   approximate leverage scores to the desired accuracy
%
%  -approxlevscores, a vector of one set of leverage score approximations


% % if testing
%out = dummy_extension(in);
%return

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
        [out.approxlevscores, out.numiters(iter)] = ...
            power_method_approx_levscores(in);
        levscoreprobs = out.approxlevscores/in.k;

        % sample according to those leverage scores
        colindices = ones(1,in.l);
        for i=1:in.l
            colindices(i) = find(cumsum(levscoreprobs) >= rand(),1);
        end
        scalingfactors = levscoreprobs(colindices).^(1/2);
        S = Id(:,colindices)*diag(1./scalingfactors);

        if in.linearkernelflag == 0
            C = in.A*S;
            [Q,~] = qr(C, 0);
        else
            C = in.A*(At*S);
            [Q,~] = qr(C, 0);
        end

        if in.linearkernelflag == 0
            W = full(Q'*in.A*Q);
        else
            Wroot = full(At*Q);
            W = Wroot'*Wroot;
        end
        
    out.timings(:, iter) = toc;
       
    % time to form Wk
    tic
        [V,D] = orderedeig(W);
        Wk = V(:, 1:in.k)*D(1:in.k, 1:in.k)*V(:,1:in.k)';
    out.timings(2, iter) = out.timings(2, iter) + toc;
    
    [out.specerr(:,iter), out.froerr(:,iter), out.trerr(:,iter)] = ...
        estnorms(in, Q, W, Wk);

end

end

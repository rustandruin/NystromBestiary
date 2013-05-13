function out = speclev_randEig2(in)
% out = speclev_randEig2(in)
%
% in is a structure with (at least) the following fields:
%  -A, an n-by-n PSD matrix, and
% -linearkernelflag, 1 if A=XX^T and X was passed in A's position, 0
%  otherwise, and
%  -k, the target rank, and 
%  -l < n, the number of columns to sample and
%  -q >= 1 the number of times to repeat the experiment and
%  -chunk, how often to reorthogonalize during leverage score approximation
%  -vareps, a precision factor in (0,1)
%
% out is a structure with the following fields:
%  - specerr, froerr, trerr: each matrices with two rows representing 
%  the spectral, frobenius, 
%  and trace norms of the errors in using q realizations of a randEig 
%  Nystrom extension of A using l columns sampled according to
%  approximate leverage score probabilities calculated using Algorithm 4 of 
%  Mahoney et al. The first row corresponds to 
%  randEig Nystrom extesions where the rank was not fixed, the second to
%  extensions where the rank was fixed
% 
%  -timings, a two row matrix of the time it took to run each experiment 
% (i.e. form Q, W, Wk), including the time to approximate the leverage scores
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

Id = speye(n);

if in.linearkernelflag
    At = in.A';
end

for iter=1:in.q
    tic
        % call whatever algorithm we're using to compute the approximate
        % leverage scores
        approxlevscores = spectral_sketch_levscores(in);
        levscoreprobs = approxlevscores/in.k;

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
            C = in.A*Q;
        else
            C = in.A*(At*S);
            [Q,~] = qr(C, 0);
            C = in.A*(At*Q);
        end

        if in.linearkernelflag == 0
            W = full(Q'*in.A*Q);
        else
            Wroot = full(At*Q);
            W = Wroot'*Wroot;
        end
        
    out.timings(:,iter) = toc;
    
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
    
    [out.specerr(:,iter), out.froerr(:,iter), out.trerr(:,iter)] = ...
        estnorms(in, C, Winv, Wkinv);
end

end

function out = qrlevscore_Nystrom(in)
% out = qrlevscore_fixedrank_Nystrom(in)
%
% in is a structure with (at least) the following fields:
%  -A, an n-by-m matrix, and
%  -k, the target rank, and 
%  -l < n, the number of columns to sample and
%  -qrlevscoreprobs, a 1-by-n vector containing the leverage score probabilities 
%  of A*A' computed using a QR factorization
%  -qrlevscorecomputationtime, the time it took to compute the leverage
%  score probabilities
%  -q >= 1 the number of times to repeat the experiment.
%
% out is a structure with the following fields:
%  -specerr, froerr, trerr: each matrices with two rows representing 
%  the spectral, frobenius, 
%  and trace norms of the errors in using q realization of the leverage
%  score based Nystrom extension to A*A'. The first row corresponds to 
%  Nystrom extensions where the rank was not fixed, the second to Nystrom
%  extensions where the rank was fixed
% 
%  -timings, a two row matrix of the time it took to run each experiment 
% (i.e. form C, Winv, Wkinv), including the time to approximate the leverage scores

% % if testing
%out = dummy_extension(in);
%return;

n = size(in.A,1);
out.specerr = zeros(2,in.q);
out.froerr = zeros(2,in.q);
out.trerr = zeros(2,in.q);
out.timings = zeros(2,in.q) + in.qrlevscorecomputationtime;

Id = speye(n);

if in.linearkernelflag % this should always be 1, for this function
    At = in.A';
end

for iter=1:in.q
    tic
        colindices = ones(1,in.l);
        for i=1:in.l
            colindices(i) = find(cumsum(in.qrlevscoreprobs) >= rand(),1);
        end
        scalingfactors = in.qrlevscoreprobs(colindices).^(1/2);
        S = Id(:,colindices)*diag(1./scalingfactors);

        if in.linearkernelflag == 0
            C = in.A*S;
        else
            C = in.A*(At*S);
        end

        W = full(S'*C);
    out.timings(:, iter) = out.timings(:, iter) + toc;
    
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
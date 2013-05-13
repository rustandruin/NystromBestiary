function out = levscore_randEig(in)
% out = levscore_fixedrank_randEig(in)
%
% in is a structure with (at least) the following fields:
%  -A, an n-by-n PSD matrix, and
% -linearkernelflag, 1 if A=XX^T and X was passed in A's position, 0
%  otherwise, and
% -k, the target rank, and 
%  -l < n, the number of columns to sample and
%  -levscoreprobs, a 1-by-n vector containing the leverage score probabilities 
%  of some dominant eigenspace of A (or approximations thereof), and
%  -levscorecomputationtime, the time it took to compute the leverage
%  score probabilities
%  -q >= 1 the number of times to repeat the experiment.
%
% out is a structure with the following fields:
%  -specerr, froerr, trerr: each matrices with two rows representing 
%  the spectral, frobenius, 
%  and trace norms of the errors in using q realization of the leverage
%  score based randEig sketch of A. The first row corresponds to 
%  randEig sketches where the rank was not fixed, the second to randEig
%  sketches where the rank was fixed
% 
%  -timings, a two row matrix of the time it took to run each experiment 
% (i.e. form Q, W, Wk), including the time to approximate the leverage scores


% % if testing
%out = dummy_extension(in);
%return

n = size(in.A,1);
out.specerr = zeros(2,in.q);
out.froerr = zeros(2,in.q);
out.trerr = zeros(2,in.q);
out.timings = zeros(2,in.q) + in.levscorecomputationtime;

Id = speye(n);

if in.linearkernelflag
    At = in.A';
end

for iter=1:in.q
    tic
        colindices = ones(1,in.l);
        for i=1:in.l
            colindices(i) = find(cumsum(in.levscoreprobs) >= rand(),1);
        end
        scalingfactors = in.levscoreprobs(colindices).^(1/2);
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
        
    out.timings(:, iter) = out.timings(:, iter) + toc;
    
    % time to form Wk
    tic
        [V,D] = orderedeig(W);
        Wk = V(:, 1:in.k)*D(1:in.k, 1:in.k)*V(:,1:in.k)';
    out.timings(2, iter) = out.timings(2, iter) + toc;
    
    [out.specerr(:,iter), out.froerr(:,iter), out.trerr(:,iter)] = ...
        estnorms(in, Q, W, Wk);
end

end

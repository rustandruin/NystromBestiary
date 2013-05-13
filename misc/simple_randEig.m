function out = simple_randEig(in)
% out = simple_randEig(in)
%
% in is a structure with (at least) the following fields:
%  -A, an n-by-n PSD matrix, and
% -linearkernelflag, 1 if A=XX^T and X was passed in A's position, 0
%  otherwise, and
%  -k, the desired rank of the approximation
%  -l < n, the number of columns to sample and
%  -q >= 1 the number of times to repeat the experiment.
%
% out is a structure with the following fields:
%  - specerr, froerr, trerr: each matrices with two rows representing 
%  the spectral, frobenius, 
%  and trace norms of the errors in using q realizations of the simple 
%  column-based (sampled uniformly at random w/o replacement)
%  randEig sketch of A, using l columns. The first row corresponds to 
%  sketches where the rank was not fixed, the second to sketches
%  where the rank was fixed
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
out.timings = zeros(2,in.q);

if in.linearkernelflag 
    At = in.A';
end

for iter=1:in.q
    tic
        p = randperm(n);
        p = p(1:in.l);
        
        if in.linearkernelflag == 0
            C = in.A(:, p);
            [Q,~] = qr(C, 0);
        else
            C = in.A*At(:, p);
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



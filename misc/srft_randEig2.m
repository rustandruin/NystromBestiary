function out = srft_randEig2(in)
% out = srft_randEig2(in)
%
% in is a structure with (at least) the following fields:
%  -A, an n-by-n PSD matrix, and
% -linearkernelflag, 1 if A=XX^T and X was passed in A's position, 0
%  otherwise, and
% -k, the target rank, and
%  -l < n, the number of columns to sample and
%  -q >= 1 the number of times to repeat the experiment.
%
% out is a structure with the following fields:
%  - specerr, froerr, trerr: each matrices with two rows representing 
%  the spectral, frobenius, 
%  and trace norms of the errors in using q realizations of the SRFT-based 
%  randEig Nystrom extension of A, using l columns. The first row corresponds to 
%  extensions where the rank was not fixed, the second to extensions
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
        [C,~] = srft(in); % for accurate timing, need to modify srft to only return C when don't need W
        [Q,~] = qr(C, 0);
        
        if in.linearkernelflag == 0
            C = in.A*Q;
            W = full(Q'*in.A*Q);
        else
            C = in.A*(At*Q);
            Wroot = full(At*Q);
            W = Wroot'*Wroot;
        end
        
    out.timings(iter) = toc;
    
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


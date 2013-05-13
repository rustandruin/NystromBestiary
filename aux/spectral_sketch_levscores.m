function approxlevscores = spectral_sketch_levscores(in)

% approxlevscores = spectral_sketch_levscores(in)
%
% approximates the leverage scores of A, filtered through rank k, 
% using Algorithm 4 of Mahoney et al.
%
% in is a structure with (at least) the following fields:
% - A, an SPSD matrix
% - linearkernelflag, chunk, k, where
%
% linearkernelflag is 0 if A = in.A is the PSD matrix, and
% 1 if the actual PSD matrix to be used is A = in.A*in.A^T
%
% chunk is the number of iterations after which to reorthogonalize 
%
% returns the approximate leverage scores and the number of iterations
% needed
%
% Assumes k <= rank(A) and A is SPSD (so A' == A)

n = size(in.A,1);
S = randn(n, 2*in.k);
% Estimate 2*log(1 + eps/10) - 1/2 from the bound in the paper as 1
% q = ceil( log(1 + sqrt(k/(k-1)) + exp(1)*sqrt(2*(n-k)/k)) );
q = 4; % if you use the above expression, for n=10000, k=30, q would be about 5

At = in.A';

if in.linearkernelflag == 0 % A is PSD
    Y = in.A*S;

    iter = 0; % keep track of how many iterations of the power method were used

    while iter < q
        iter = iter + 1;
        % reorthogonalize every chunk steps
        if (rem(iter,in.chunk)==0)
            [Q,R] = qr(Y, 0);
            Y = Q; % With probability 1, Q still has l columns
        end
        Y = in.A*(in.A*Y); % since A is PSD, Y = A*A^T*Y
    end
else % A = formalarg(A)*formalarg(A)^T
     Y = in.A*(At*S);

    iter = 0; % keep track of how many iterations of the power method were used

    while iter < q
        iter = iter + 1;
        % reorthogonalize every chunk steps
        if (rem(iter,in.chunk)==0)
            [Q,R] = qr(Y, 0);
            Y = Q; % With probability 1, Q still has k columns 
        end
        Y = in.A*(At*Y); 
        Y = in.A*(At*Y);
    end
end

[U,~,~] = svds(Y, in.k);
approxlevscores = sum(U.^2, 2)';
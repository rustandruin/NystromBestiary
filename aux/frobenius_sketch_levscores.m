function approxlevscores = frobenius_sketch_levscores(in)
% approxlevscores = frobenius_sketch_scores(in)
%
% returns additive approximations to the rank-k levscores of the columns
% of A. computed  using Algorithm 5 from Mahoney et al. "Fast 
% approximation of matrix coherence and statistical leverage"
% http://arxiv.org/abs/1109.3843
%
% in is a structure with (at least) the following fields:
% - A, an SPSD matrix
% - chunk, k, where
%
% chunk is the number of iterations after which to reorthogonalize 
%
% Assumes k <= rank(A) and A is SPSD (so A' == A)


n = size(in.A,2); % right dimension for both A PSD and A = Aformal*Aformal^T
%r = k + ceil(10*k + 1); % the algorithm suggests r = k + ceil(10*k/vareps + 1)
r = 2*in.k; 
r = min(n,r); % Shouldn't be necessary, since hopefully k << n/2
% could use SRFT samples here instead of Gaussians
S = randn(n, r);
B = in.A*S;

[Q,~] = qr(B, 0);
[U,~,~] = svds(Q'*in.A, in.k);
approxlevscores = sum((Q*U).^2, 2)';

end


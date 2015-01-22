function [C,W] = srft(in)
%
% [C, W] = srft(in)
%
% Given a PSD matrix A of size n, k<=n, and l<=n,
% returns  C = A*S and W = S'*A*S, where S is an n*l (real) SRFT matrix, 
% that is, S = sqrt(n/l)*D F^t R^t, where D is a diagonal matrix of random 
% signs, F is the real FFT matrix, and R restricts from n coordinates to l
%
% in is a structure with (at least) the following fields:
%  -A, an n-by-n PSD matrix, and
% -linearkernelflag, 1 if A=XX^T and X was passed in A's position, 0
%  otherwise, and
% -k, the target rank, and
%  -l < n, the number of columns to sample
%

n = size(in.A,1);
colindices = randperm(n);
colindices = colindices(1:in.l);
d = 1-2*(rand(1,n) > .5);

% compute A*D*F' "efficiently" using realfft, the fact that A is symmetric, 
% and that A*D is 10 or so times slower than using bsxfun
% my implementation of the realfft is about 5 times slower than the regular
% fft

if in.linearkernelflag == 0 % A is PSD
    Y = realfft(full(bsxfun(@times,in.A,d)'));
    C = Y(colindices, :)';
else % A = formalarg(A)*formalarg(A)^T
    % assume that formalarg(A) fits in memory when made full
    Y = realfft(full(bsxfun(@times, in.A, d')))';
    C = in.A*Y(:, colindices);
end

Y = realfft(full(bsxfun(@times,C,d'))); % compute F*D*C efficiently 
W = Y(colindices, :);

end

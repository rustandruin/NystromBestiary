function [approxlevscores, iter] = power_method_approx_levscores(in)
% [approxlevscores, iter] = power_method_approx_levscores(in)
%
% approximates the leverage scores of A, filtered through rank k, 
% using the power method. The number of iterations is determined by testing
% that the approximations don't change by tol, or that a maximum number of
% iterations has been reached.
%
% in is a structure with (at least) the following fields:
% - A, an SPSD matrix
% - linearkernelflag, chunk, k, tol, maxiters, where
%
% linearkernelflag is 0 if A = in.A is the PSD matrix, and
% 1 if the actual PSD matrix to be used is A = in.A*in.A^T
%
% chunk is the number of iterations after which to reorthogonalize and to check
% the leverage scores for convergence
%
% tol determines the tolerance for convergence: when the inf norm distance
% of the approximate leverage scores is smaller than this, convergence 
% has been achieved
%
% maxiters is the maximum number of iterations to use
%
% returns the approximate leverage scores and the number of iterations
% needed
%
% Assumes k <= rank(A) and A is SPSD (so A' == A)

n = size(in.A,1);
l = in.k;
S = randn(n, l);

if in.linearkernelflag == 0 % A is SPSD
    Y = in.A*S;

    iter = 0; % keep track of how many iterations of the power method were used

    approxlevscores = zeros(1,n);

    while iter < in.maxiters
        iter = iter + 1;
        if (rem(iter,in.chunk)==0)
            [Q,R] = qr(Y, 0);
            Y = Q; % With probability 1, Q still has k columns 
            oldlevscores = approxlevscores;
            approxlevscores = sum(Q.^2,2)';
            if ( norm(approxlevscores - oldlevscores, Inf) < in.tol)
                return;
            end
        end
        Y = in.A*(in.A*Y); % since A is PSD, Y = A*A^T*Y
    end
else % A = formalarg(A)*formalarg(A)^T
     Y = in.A*(in.A'*S);

    iter = 0; % keep track of how many iterations of the power method were used

    approxlevscores = zeros(1,n);

    while iter < in.maxiters
        iter = iter + 1;
        if (rem(iter,in.chunk)==0)
            [Q,R] = qr(Y, 0);
            Y = Q; % With probability 1, Q still has k columns 
            oldlevscores = approxlevscores;
            approxlevscores = sum(Q.^2,2)';
            if ( norm(approxlevscores - oldlevscores, Inf) < in.tol)
                return;
            end
        end
        Y = in.A*(in.A'*Y);
        Y = in.A*(in.A'*Y);
    end
end

if (iter == in.maxiters)
    fprintf('Terminated with current approximate leverage scores after %i iterations\n', in.maxiters);
end

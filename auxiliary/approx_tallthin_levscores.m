function [proj1time, proj2time, overalltime, levscores] = approx_tallthin_levscores(A, vareps, method)
% levscores = approx_tallthin_levscores(A, vareps, method)
%
% Takes a full-rank tall thin matrix A in R^{n times d} and an error parameter
% vareps in (0, 1/2]
%
% Returns approximations to the leverage scores of the rows of A, using
% Algorithm 1 in Mahoney et al.
[n, d] = size(A);

delta = 1/10; % failure probability of the FJLT
% the formula for r1 is taken from Gittens-Boutsidis, Improved Matrix ...
r1 = ceil(8/3*vareps^(-2)*log(3/delta*d)*(sqrt(d) + sqrt(8*log(3/delta*n)))^2);
r1 = ceil(vareps^(-2)*log(1/delta*d)*(sqrt(d) + sqrt(log(1/delta*n)))^2);
% the formula for r2 is taken from Mahoney et al.
r2 = ceil(vareps^(-2)*(12*log(n) + 6*log(1/delta))); 
r2 = ceil(vareps^(-2)*(log(n) + log(1/delta))); 

fprintf('r1/n = %.3f, r2/d = %.3f\n, r1 = %d, r2 = %d\n', r1/n, r2/d, r1, r2);

if(r1 > n)
    fprintf('Warning: r1 > n\n');
end
if (r2 > d)
    fprintf('Warning: r2 > d\n');
end


tic
signs = 1-2*(rand(n,1) > .5);
rowindices = randperm(n);
rowindices = rowindices(1:min(n,r1));
Y = realfft(full(bsxfun(@times,signs,A)));
Y = Y(rowindices, :);
proj1time = toc;

overalltime = proj1time;
tic;
switch method
 case 'svd'
    [~, SY, VY] = svd(Y, 'econ');
    Yp = (A*VY)/SY;
 case 'qr'
    [Q,R] = qr(Y, 0);
    Yp = A/R;
end
overalltime = overalltime + toc;

tic
randnumbers = rand(1,r2*d);
poslocs = randnumbers < 1/6;
neglocs = randnumbers >= 1/6 & randnumbers < 1/3;

JLT = spalloc(r2, d, length(poslocs)+length(neglocs));
JLT(poslocs) = sqrt(3/r2);
JLT(neglocs) = -sqrt(3/r2);
proj2time = toc;

tic
Omega = Yp*JLT';
overalltime = overalltime + toc + proj2time;

tic
levscores = sum(Omega .* Omega, 2);
% normalize them to sum to d
levscores = d*levscores/sum(levscores);
overalltime = overalltime + toc;

end

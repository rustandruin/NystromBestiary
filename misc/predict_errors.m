function bounds = predict_errors(saveddata, delta)

% Computes the error guarantees for several different bounds in the
% literature, in the form of relative error bounds for approximating
% a kernel matrix using l samples, relative to the actually observed
% rank k errors. 
%
% To use some of the bounds, would need much more than the l samples that 
% were actually used if we applied the bounds without change, so 
% ignore the constant factors in 
% the number of samples needed, and where coherence occurs in the 
% expressions for l, assume it is 1. The parameter epsilon in error bounds
% is determined by chosing the value of epsilon which when substituted into
% the expression for l yields the actual value of l used.
%
% predict_errors(savedata, delta)
% - saveddata, a structure containing the results of a run of Nystrom
%   experiments, in the form returned by generate_dataset
% - delta, a failure probability
%
% Outputs:
% -bounds contains the resulting relative errors and the computed effective
%  epsilons (so can check that the effective epsilons are smaller than one
%  where necessary: if they aren't then the bounds do not apply)
%
% bounds.empirical* -- the empirically observed errors stored in saveddata
% bounds.{lev*,srft*,gauss*,nystrom*} -- the bounds from our paper
% bounds.drineas* --- the bounds from "On the Nystrom method for Approximating
%   A Gram Matrix ..." by Drineas and Mahoney
% bounds.belabbastrbound -- the trace bound from "On landmark selection..." by
%  Belabbas and Wolfe
% bounds.talwalkar* -- the bounds from "Ensemble Nystrom Method"
%   by Kumar, Mohri, Talwalkar

if saveddata.in.linearkernelflag
    A = saveddata.in.A*saveddata.in.A';
else
    A = saveddata.in.A;
end
n = size(A,1);
l = saveddata.in.lvals;
k = saveddata.in.k;

optspecerr = saveddata.optspecerr;
optfroerr = saveddata.optfroerr;
opttrerr = saveddata.opttrerr;
sumsquareddiagonals = sum(diag(A).^2);

% Empirical results
bounds.empiricalnystrom2 = mean(saveddata.simpleData.specerr(1,:));
bounds.empiricalnystromF = mean(saveddata.simpleData.froerr(1,:));
bounds.empiricalnystromtr = mean(saveddata.simpleData.trerr(1,:));

bounds.empiricalsrft2 = mean(saveddata.srftData.specerr(1,:));
bounds.empiricalsrftF = mean(saveddata.srftData.froerr(1,:));
bounds.empiricalsrfttr = mean(saveddata.srftData.trerr(1,:));

bounds.empiricallev2 = mean(saveddata.levscoreData.specerr(1,:));
bounds.empiricallevF = mean(saveddata.levscoreData.froerr(1,:));
bounds.empiricallevtr = mean(saveddata.levscoreData.trerr(1,:));

bounds.empiricalgauss2 = mean(saveddata.gaussianData.specerr(1,:));
bounds.empiricalgaussF = mean(saveddata.gaussianData.froerr(1,:));
bounds.empiricalgausstr = mean(saveddata.gaussianData.trerr(1,:));

% Drineas et al. bound 
eta = 1 + sqrt(log(1/delta));
epsilon = (k*eta^2/l)^(1/4);
bounds.drineas2bound = (optspecerr + epsilon*sumsquareddiagonals)/...
    bounds.empiricalnystrom2;
bounds.drineasFbound = (optfroerr + epsilon*sumsquareddiagonals)/...
    bounds.empiricalnystromF;

% Belabbas and Wolfe bound
bounds.belabbastrbound = 1/delta*(n-l)/n*trace(A)/bounds.empiricalnystromtr;

% Matrix coherence bound
r = (norm(A, 'fro')/svds(A,1))^2; % an underestimate of the rank of A
bounds.matcohrankestimate = r; % exact recovery holds if r > l

% Talwalkar et al. (sampling methods paper) bounds
% - Kmax, the maximum entry on the diagonal of the input matrix A
% - dmax, the maximum entry in the distance matrix associated with A:
%   max of sqrt(A_{ii} + A_{jj} -2A_{ij}) over i,j
diagentries = diag(A);
auxmat = bsxfun( @plus, diagentries, diagentries');
auxmat = sqrt(auxmat - 2*A);
Kmax = max(diagentries);
dmax = max(auxmat(:));

beta = 1-1/(2*max(l, n-l));
prod1 = (n-l)/(n-1/2)*1/beta*log(1/delta);
bounds.talwalkar2bound = (optspecerr + ...
                  n/sqrt(l)*Kmax*(1 + sqrt(prod1)*dmax/Kmax^(1/2)))/...
                  bounds.empiricalnystrom2;
bounds.talwalkarFbound = (optfroerr + ...
              (k/l)^(1/4)*n*Kmax*sqrt(1 + sqrt(prod1)*dmax/Kmax^(1/2)))/...
              bounds.empiricalnystromF;

% The following bounds only apply if the corresponding epsilon is smaller
% than 1

% Lemma 1: leverage
epsilon = sqrt(k*log(k/delta)/l); % doesn't account for coherence
bounds.lev2bound = (optspecerr + epsilon^2*opttrerr)/bounds.empiricallev2;
bounds.levFbound = (optfroerr + epsilon^2*opttrerr)/bounds.empiricallevF;
bounds.levtrbound = (1 + epsilon^2)*opttrerr/bounds.empiricallevtr;
bounds.levepsilon = epsilon; 

% Lemma 2: fourier
epsilon = (sqrt(k) + sqrt(log(n/delta)))^2*log(k/delta)/l;
bounds.srft2bound = ((1 + 1/(1 - sqrt(epsilon))*(1 + log(n/delta)^2/l))*optspecerr + ...
             log(n/delta)/((1-sqrt(epsilon))*l)*opttrerr)/...
             bounds.empiricalsrft2;
bounds.srftFbound = (optfroerr + epsilon*opttrerr)/...
                    bounds.empiricalsrftF;
bounds.srfttrbound = (1 + epsilon)*opttrerr/bounds.empiricalsrfttr;
bounds.srftepsilon = epsilon;

% Lemma 3: gaussian (using linear number of samples)
epsilon = sqrt((l/k)-1)^(-1/2);
bounds.gauss2bound = ((1 + epsilon^2)*optspecerr + ...
    epsilon^2/k*opttrerr)/bounds.empiricalgauss2;
bounds.gaussFbound = (optfroerr + ...
     epsilon/sqrt(k)*opttrerr + epsilon^2*sqrt(optspecerr*opttrerr) + ...
     epsilon^2*sqrt(log(k)/k)*optspecerr)/...
     bounds.empiricalgaussF;
bounds.gausstrbound = (1 + epsilon^2)*opttrerr/...
    bounds.empiricalgausstr;
bounds.gaussepsilon = epsilon;

% Lemma 4: nystrom
epsilon = sqrt(k*log(k/delta)/l); % doesn't account for coherence
bounds.nystrom2bound = (1 + n/((1-epsilon)*l))*optspecerr/bounds.empiricalnystrom2;
%bounds.nystrom2bound = (optspecerr + 1/(delta^2*(1-epsilon))*opttrerr)/bounds.empiricalnystrom2;
bounds.nystromFbound = (optfroerr + ...
    1/(delta*sqrt(1-epsilon))*opttrerr)/bounds.empiricalnystromF;
bounds.nystromtrbound = (1 + 1/(delta^2*(1-epsilon)))*opttrerr/bounds.empiricalnystromtr;
bounds.nystromepsilon = epsilon;

end

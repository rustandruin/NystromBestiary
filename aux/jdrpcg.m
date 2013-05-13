function varargout=jdrpcg(varargin)
%JDRPCG computes eigenpairs of a square real symmetric matrix or operator.
% 
%  Lambda = JDRPCG(A)
%  returns the K smallest eigenvalues in a K vector Lambda. Here 
%  K=min(5,N) (unless K has been specified), where N=size(A,1). JDRPCG(A) 
%  (without output argument) displays the K eigenvalues.
%
%  [X,Lambda] = JDRPCG(A) 
%  returns the eigenvectors in the N by K matrix X and the eigenvalues
%  in the K by K diagonal matrix Lambda: A*X = X*Lambda. 
%
%  ... = JDRPCG('Afun')
%  ... = JDRPCG('Afun',N)
%  The first input argument is either a square matrix (which can be full
%  or sparse, symmetric or nonsymmetric, real or complex), or a string
%  containing the name of an M-file which applies a linear operator to a
%  given column vector. In the latter case, the M-file must return the 
%  order N of the problem with N = Afun([],'dimension') or N must be
%  specified in the list of input arguments. For example, EIGS('fft',...)
%  is much faster than EIGS(F,...) where F is the explicit FFT matrix.
%
%
%  Other input arguments are optional and can be given in practically any
%  order:
%  ... = JDRPCG(...,K,SIGMA,OPTIONS,M) 
%  where
%
%      K         An integer, the number of eigenvalues desired.
%      SIGMA     A scalar shift.
%      OPTIONS   A structure containing additional parameters.
%      M         A string or a matrix that specifies the preconditioner.
%
%  If K is not specified, then K = MIN(N,5) eigenvalues are computed.
%
%  If SIGMA is specified, it is used as a lower bound 
%  on the smallest eigenvalue, whose knowledge improves algorithm
%  efficiency when the distance from this lower bound to the samllest
%  eigenvalue is not much larger than the distance between the two 
%  smallest eigenvalues. 
%
%  The OPTIONS structure specifies certain parameters in the algorithm.
%
%    Field name         Parameter                             Default
%
%    OPTIONS.Tol        Convergence tolerance for              1e-8 
%                         the residual norm of each eigenpair
%
%    OPTIONS.jmin       Minimum dimension search subspace.     7
%    OPTIONS.jmax       Maximum dimension search subspace.     14
%
%    OPTIONS.MaxIt      Maximum number of matrix vector        5000
%                                            product with A. 
%    OPTIONS.v0         Starting vector.                       ones
%
%    OPTIONS.Disp       If 'Disp'=1 then,                      0
%                           the input parameters, 
%                           the residual size at each step,
%                           and the eigenvalues at detection 
%                         are displayed, and      
%                         the convergence history is plotted.
%
%    OPTIONS.Precond    Preconditioner.                        LU=[[],[]].
%
%  For instance
%
%    OPTIONS=STRUCT('Tol',1.0e-10,'Precond',M);
%
%  changes the convergence tolerance to 1.0e-10 and takes the preconditioner 
%    defined in M.m if M is the string 'M', 
%    or M = L*U if M is an n by 2*n matrix: M = [L,U].
%
%  The preconditoner can be specified in the OPTIONS structure, but also
%  in the argument list:
%   ... = JDRPCG(...,K,SIGMA,M,OPTIONS) 
%   ... = JDRPCG(...,K,SIGMA,L,U,OPTIONS) 
%   ... = JDRPCG(...,K,SIGMA,'M',OPTIONS)
%   ... = JDRPCG(...,K,SIGMA,'L','U',OPTIONS)
%  as an N by N matrix M (then M is the preconditioner), or an N by 2*N 
%  matrix M (then  L*U is the preconditioner, where  M = [L,U]), or as
%  N by N matrices L and U (then  L*U is the preconditioner), or as one 
%  or two strings containing the name of  M-files ('M', or 'L' and 'U')
%  which apply a linear operator to a given column vector.
%
%
%  [X,Lambda,HISTORY] = JDRPCG(A,...) 
%  returns also the convergence history. HISTORY is an array of 5 columns: 
%     HISTORY(:,1) is the index of the eigenpair being searched;
%     HISTORY(:,2) is the number of JD iterations already completed;
%     HISTORY(:,3) is the cumulative number of multiplications by A;
%     HISTORY(:,4) is the norm of the residual;
%     HISTORY(:,5) is the current approximation of the searched eigenvalue.
%  Here, a JD iteration is an iteration of the main outer loop of the 
%  algorithm, and the number of such iterations is equal to the number of
%  call to the PCG inner solver, see the technical report referenced below.
%
%  JDRPCG(without input arguments) lists the options and the defaults.
%
%---------------------------------------------------------------------------
%   Modified version of the JDQR package by Gerard Sleijpen.
%   Copyright (c) Gerard Sleijpen 1998 & Yvan Notay 2000, 2005.
%---------------------------------------------------------------------------
%   The algorithm used is Jacobi-Davidson with regular preconditioning and
%   CG inner iterations, as described in
%
%     Y. Notay,
%     Inner iterations in eigenvalue solvers,
%     Report GANMN 05-01, ULB, Brussels, Belgium, 2005.
%
%   Some implementations features from the following report are also used
%
%     Y. Notay, 
%     Combination of Jacobi-Davidson and conjugate gradients for the partial
%     symmetric eigenproblem, Report GANMN 00-01, ULB, Brussels, Belgium, 2000.
%     (appeared in Numer. Lin. Alg. Appl., vol. 9, pp. 21-44, 2002).
%
%
% Both reports are available at  http://homepages.ulb.ac.be/~ynotay/
%
% Other references and links about the Jacobi-Davidson method can be found at
% http://www.math.uu.nl/people/sleijpen/ 
% http://www.math.ruu.nl/people/vorst/
%----------------------------------------------------------------------------

global Qschur Lambda nm_operations history

if nargin==0, ShowOptions, return, end

%%% Read/set parameters
[n,nselect,sigma,...
   jmin,jmax,tol,maxit,V,SHOW] = ReadOptions(varargin{1:nargin}); 

%%% Initiate global variables
Qschur = zeros(n,0); Lambda = []; rr=[];
nm_operations = 0; history = []; histjd=[];

%%% Return if eigenvalueproblem is trivial
if n<2
  if n==1, Qschur=1; Lambda=MV(1); end
  if nargout == 0, eigenvalue=Lambda, else
  [varargout{1:nargout}]=output(history); end, 
return, end

String = '\r#it=%i #MV=%i dim(V)=%i |r_%i|=%6.1e ';
time = clock;

%%% Initialize V, W:
W=MV(V); M=W'*V;
j=1; k=0;
nit=0; del=-1; FIG=get(0,'CurrentFigure'); 
nr=Inf;
if isreal(sigma)&&isnumeric(sigma)&&~isempty(sigma), eta0=sigma(1); 
else eta0=Inf; end

%%% The JD loop (Standard)
%%%    V orthogonal, V orthogonal to Qschur
%%%    V*V=eye(j), Qschur'*V=0, 
%%%    W=A*V, M=V'*W
%%%
while (k<nselect) && (nm_operations <= maxit) 

   %%% Compute approximate eigenpair and residual
 [UR,S]=SortSchur(M);

   %%% Select approximate eigenpair from subspace if a converged eigenvector
   %%% has not been found during CG  
 if nr >= tol
   y=UR(:,1); theta=S(1,1); u=V*y; w=W*y; r=w-theta*u;  nr=norm(r);
   eta0=min(eta0,theta);
   % explicit check of residual norm if convergence reached
     if nr < tol
       w=MV(u); theta=w'*u;  r=w-theta*u;  nr=norm(r);
     end
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   history=[history;k+1,nit,nm_operations,nr,theta];               %%%
   if SHOW, fprintf(String,nit,nm_operations,j,nit,nr)             %%%
   histjd=[histjd;nr,nm_operations];  end                          %%%
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 end

 %%% Proceed converged eigenpairs
 while nr < tol

      %%% Expand the partial Schur form
      Qschur=[Qschur,u]; 
      Lambda=[Lambda,theta]; rr=[rr,nr];  k=k+1;
      if SHOW, ShowLambda(theta,k), end

      if k>=nselect, break, end

      if j==1
          V = ones(n,1)+0.1*rand(n,1); V=RepGS(Qschur,V);         
          W=MV(V); M=W'*V;
      else
         j=min(jmin+1,j-1);
         J=[2:j+1]; UR=UR(:,J); 
         M=S(J,J); V=V*UR; W=W*UR; 
      end

        % already compute next residual
        u=V(:,1); w=W(:,1); theta=M(1,1); r=w-theta*u;  nr=norm(r);
        % explicit check of residual norm if convergence reached
        if nr < tol
           w=MV(u); theta=w'*u;  r=w-theta*u;  nr=norm(r);
        end
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   history=[history;k+1,nit,nm_operations,nr,theta];               %%%
   if SHOW, fprintf(String,nit,nm_operations,j,nit,nr)             %%%
   histjd=[histjd;nr,nm_operations];  end                          %%%
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

     del=-1;
     eta0 = max([eta0,Lambda(1:k)]);

   UR=eye(j); S=M; y=[1 ; zeros(j-1,1)];

 end % nr<tol

if k>=nselect || nm_operations==maxit, break, end

 %%% Check for shrinking the search subspace
 if j>=jmax
      j=jmin; J=[1:j]; UR=UR(:,J);
      M=S(J,J); V=V*UR; W=W*UR; y=[1 ; zeros(j-1,1)];
 end % if j>=jmax
 
 %%% Solve correction equation
             eta=eta0;
             if (j>1)
                delo=del;
                del=S(2,2)-S(1,1);
                if  ( nr<del && abs(del/delo-1)<0.1 ), eta=theta; end
             end

if k>0, r=RepGS(Qschur,r,0); end

[v,rk,iet]=Solve_pce(eta,u,r,tol,maxit-nm_operations-1,theta,nr);
      nit=nit+1;

      if rk < tol   % convergence has been reached during CG
         u=u-v;
         if k>0, u=RepGS(Qschur,u,0); end
         u=u/norm(u); 
         w=MV(u); theta=w'*u;  r=w-theta*u;  nr=norm(r);
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   history=[history;k+1,nit,nm_operations,nr,theta];               %%%
   if SHOW, fprintf(String,nit,nm_operations,j,nit,nr)             %%%
   histjd=[histjd;nr,nm_operations];  end                          %%%
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      end

      if iet==3 && eta<=eta0 && k>0          % eta0 too large
         eta0=min([eta0,Lambda(1:k)]); end

 %%% Expand the subspaces of the interaction matrix  
 if nr >= tol || k < nselect-1
        [v,yy]=RepGS([Qschur,V],v);
           w=MV(v);
           z= V'*w;
           M=[M , z ; z', v'*w]; 
           V=[V,v]; W=[W,w]; j=j+1;
        if nr < tol
           y=[y;0]-yy(k+1:k+j);
           y=y/norm(y);
           z=M*y-theta*y;
           M=M-z*y'-y*z';
        end
 end

end % while


time_needed=etime(clock,time);

%-------------- display results ----------------------------
if SHOW && size(history,1)>0
   if isempty(FIG), FIG=1; end, figure(FIG)

   StringT=sprintf('Jacobi-Davidson with jmin=%g, jmax=%g, residual tolerance %g.',...
           jmin,jmax,tol); 
   StringX=sprintf('Correction equation solved with PCG');   
   date=fix(clock);
   String=sprintf('\n%2i-%2i-%2i, %2i:%2i:%2i',date(3:-1:1),date(4:6));

     t=history(:,3);
     plot(histjd(:,2),log10(histjd(:,1)),'b*',...
          t,log10(history(:,4)),'b-',t,log10(tol)+0*t,'k:' )
     legend('log_{10} || r_{#MV} ||_2 at JD steps','log_{10} || r_{#MV} ||_2 during PCG'), 
     title(StringT)
     xlabel([StringX,String])
   drawnow

   str1=num2str(abs(k-nselect)); str='s';
   if k<nselect,
     if k==0, str1='any'; str=''; elseif k==nselect-1, str1='one'; str=''; end
     fprintf('\n\nFailed detection of %s eigenpair%s.',str1,str)
   end

   if k>0, 
      Str=[];
      for j=1:size(Lambda,2);
          Str=[Str,sprintf('\n%1s','')]; 
          Str=[Str,sprintf(' %+22.15e',Lambda(j))]; 
          Str=[Str,sprintf(' %30.1e',rr(j) )];
      end
    Str=[sprintf('\nDetected eigenvalues:     Corresponding residual norm:%15s',''),Str];
    fprintf('%s\n',Str)
   else, fprintf('\n'); 
   end

   Str='time_needed';                              DispResult(Str,eval(Str))
   if (k>0)
      I=eye(k); Str='norm(X''*X-I)';     DispResult(Str,norm(Qschur'*Qschur-I))  
   end
   fprintf('\n\n')

end

if nargout == 0, if ~SHOW, eigenvalues=Lambda, end, return, end
[varargout{1:nargout}]=output(history);

return
%========== OUTPUT =========================================================
function varargout=output(history)

global Qschur Lambda

if nargout==1,    varargout{1}=Lambda; return
else  
    varargout{1}=Qschur; varargout{2}= Lambda;
    if nargout == 3, varargout{3}=history; end
end
return

%===========================================================================
%===== SOLVE CORRECTION EQUATION ===========================================
%===========================================================================
function [t,rk,iet]=Solve_pce(eta,u,r,tol,max_it,theta,g0)

global Qschur L_precond

         [t,rk,iet] = cgjd(eta,[Qschur,u],u,r,tol,max_it,theta,g0);

return
%=======================================================================
function [t,rk,iet] = cgjd(eta,Q,u,g,tol,max_it,theta,g0)
%
global nm_operations history
%
% Parameters
sg=sqrt(0.1); trs=15; pl=1; st=0.1;
%
n = size(g,1); iet=0;
bk=0; rho=1; gk=g0; rk=gk; thk=theta;

g_trs=max(sg*g0,tol); g_s=g_trs;
gkmin=g0; tmin=zeros(n,1);

if max_it==0 || tol>=1, t=g;  return, end

i=size(history,1);
nsearch=history(i,1);
nit=history(i,2);

for i = 1 : max_it
    if (i > 1)
      gk=norm(g);
      if ( gk <= gkmin )
         gkmin=gk;
         tmin=t;
      end
      if ( gk < g_trs ) 
          update=1;
          if ( gk < g_s )           
              t_est=norm(t); beta_est=abs(theta-eta+bk); update=0;
              if (gk < st*g0), g_s=tol ;else, g_s=max(st*g0,tol); end
          end
          thk=beta_est/(1+t_est^2); 
          rk=sqrt( gk^2*(1+t_est^2) + (t_est*beta_est)^2 )/(1+t_est^2);
          history=[history;nsearch,nit,nm_operations,rk,thk+eta];
          if (rk < tol && update) 
             t_est=norm(t); beta_est=abs(theta-eta+bk);
             rk=sqrt( gk^2*(1+t_est^2) + (t_est*beta_est)^2 )/(1+t_est^2);
          end
          if (rk < tol), return, end
          if ( gk < trs*beta_est*t_est/sqrt(1+t_est^2) && ...
                    beta_est*t_est/(1+t_est^2) > 0.5*tol )
                iet=1; return, end          
          if ( gk > pl*gkmin && beta_est*t_est/(1+t_est^2) > 0.5*tol )
              t=tmin; iet=2; return, end
      end
    end

   w = SolvePrecond(g);
   w=w-Q*(Q'*w);
   
   rho1 = rho;
   rho = g' * w;
   if ((rho <= 0) || isinf(rho))
      error(' The preconditioner is not positive definite')
   end
   if (i == 1)
      d = w;
   else
      beta = rho / rho1;
      d = w + beta * d;
   end

   v=MV(d)-eta*d;
   v=v-u*(u'*v);
   
   alpha=v' * d;
   if (alpha <= 0), iet=3; if i==1, t=d; else, t=tmin; end, return, end
   gam = rho / alpha;
   if (i == 1)   
     t = gam * d; tmin=t;
   else
     t = t + gam * d;
   end
   g = g - gam * v;

   bk=bk-rho*gam;      

end                                % for i = 1 : maxit
return
%======================================================================
%========== BASIC OPERATIONS ==========================================
%======================================================================
function v=MV(v)

global A_operator nm_operations

if ischar(A_operator)
  for i=1:size(v,2), v(:,i) = feval(A_operator,v(:,i)); end
else
  v = A_operator*v; 
end

nm_operations = nm_operations+size(v,2);

return
%----------------------------------------------------------------------
function y=SolvePrecond(y,flag)
% Action preconditioner

  global type_preconditioner L_precond U_precond P_precond

  if nargin==1
    switch type_preconditioner
      case 0,     
      case 1,     y=feval(L_precond,y); 
      case 2,     y=feval(L_precond,y,'preconditioner');
      case 3,     y=feval(L_precond,y); 
                  y=feval(U_precond,y);
      case 4,     y=feval(L_precond,y,'L'); 
                  y=feval(L_precond,y,'U');
      case 5,     y=L_precond\y;
      case 6,     y=U_precond\(L_precond\y);
      case 7,     y=U_precond\(L_precond\(P_precond*y));
    end
  elseif strmatch(flag,'U')
    switch type_preconditioner
      case 3,        y=feval(U_precond,y);
      case {1,2,4},  y=feval(L_precond,y,'U');
      case {6,7},    y=U_precond\y;
    end 
  elseif strmatch(flag,'L')
    switch type_preconditioner
      case 3,        y=feval(L_precond,y,'L');
      case {1,2,4},  y=feval(L_precond,y,'L');
      case 6,        y=L_precond\y;
      case 7,        y=L_precond\(P_precond*y);
    end
  end
    
return
%----------------------------------------------------------------------
function  r=OrtProj(Q,r);

   if ~isempty(Q), 
      r=r-Q*(Q'*r);
   end 

return
%----------------------------------------------------------------------
%=======================================================================
%========== Orthogonalisation ==========================================
%=======================================================================
function [v,y]=RepGS(V,v,gamma)
% [v,y]=REP_GS(V,w)
% If V orthonormal then [V,v] orthonormal and w=[V,v]*y;
% If size(V,2)=size(V,1) then w=V*y;
%
% The orthonormalisation uses repeated Gram-Schmidt
% with the Daniel-Gragg-Kaufman-Stewart (DGKS) criterion.
%
% [v,y]=REP_GS(V,w,GAMMA)
% GAMMA=1 (default) same as [v,y]=REP_GS(V,w)
% GAMMA=0, V'*v=zeros(size(V,2)) and  w = V*y+v (v is not normalized).


% coded by Gerard Sleijpen, August 28, 1998

if nargin < 3, gamma=1; end

[n,d]=size(V);

if size(v,2)==0, y=zeros(d,0); return, end

nr_o=norm(v); nr=eps*nr_o; y=zeros(d,1);
if d==0
  if gamma, v=v/nr_o; y=nr_o; else, y=zeros(0,1); end, return
end

y=V'*v; v=v-V*y; nr_n=norm(v); ort=0;

while (nr_n<0.5*nr_o && nr_n > nr)
  s=V'*v; v=v-V*s; y=y+s;
  nr_o=nr_n; nr_n=norm(v);     ort=ort+1;
end

if nr_n <= nr, if ort>2, disp(' dependence! '), end
  if gamma  % and size allows, expand with a random vector
    if d<n, v=RepGS(V,rand(n,1)); y=[y;0]; else, v=zeros(n,0); end
  else, v=0*v; end
elseif gamma, v=v/nr_n; y=[y;nr_n]; end

return
%=======================================================================
%============== Sorts Schur form =======================================
%=======================================================================
function [Q,S]=SortSchur(A)
%[Q,S]=SortSchur(A)
%  A*Q=Q*S with diag(S) in ascending order.
%  (A real symmetric)
%
  l=size(A,1);
  if l<2, Q=1;S=A; return, end

%%%------ compute schur form -------------
  [Q,S]=eig(A); %% A*Q=Q*S, Q'*Q=eye(size(A));

%%%------ find order eigenvalues ---------------
  [D,I] = sort(diag(S)); 

%%%------ reorder schur form ----------------
  S=S(I,I);
  Q=Q(:,I);

return
%=======================================================================
%======= SET PARAMETERS ================================================
%=======================================================================
function [n,nselect,sigma,...
         jmin,jmax,tol,maxit,V,SHOW] = ReadOptions(varargin)
% Read options and set defaults

global A_operator type_preconditioner L_precond U_precond P_precond

A_operator = varargin{1};  

%%% determine dimension
if ischar(A_operator)
  n=-1;
  if exist(A_operator) ~=2
    msg=sprintf('  Can not find the M-file ''%s.m''  ',A_operator);
    errordlg(msg,'MATRIX'),n=-2;
  end
  if n==-1, eval('n=feval(A_operator,[],''dimension'');','n=-1;'), end
else
  [n,n] = size(A_operator);
  if any(size(A_operator) ~= n)
    msg=sprintf('  The operator must be a square matrix or a string.  ');
    errordlg(msg,'MATRIX'),n=-3;
  end
end

%%% defaults
jmin    =  7;
jmax    = 14; 
tol     = 1e-8; 
maxit   = 5000;
V       = zeros(0,0);
SHOW    = 0;
 
options=[]; sigma=[]; varg=[]; L_precond = []; U_precond = []; P_precond = [];
for j = 2:nargin
  if isstruct(varargin{j})
    options = varargin{j};
  elseif ischar(varargin{j}) 
    if  isempty(L_precond)
      L_precond=varargin{j};
    elseif isempty(U_precond)
      U_precond=varargin{j};
    end
  elseif length(varargin{j}) == 1
    varg = [varg,varargin{j}];
  elseif min(size(varargin{j}))==1 
    sigma = real(varargin{j}); if size(sigma,1)==1, sigma=sigma'; end 
  elseif isempty(L_precond)
    L_precond=varargin{j};
  elseif isempty(U_precond)
    U_precond=varargin{j};
  elseif isempty(P_precond)
    P_precond=varargin{j};
  end
end


[s,I]=sort(varg); I=flipdim(I,2);  
J=[]; j=0; 
while j<length(varg)
  j=j+1; jj=I(j); s=varg(jj);
  if isreal(s) && (s == fix(s)) && (s > 0)
    if n==-1
      n=s; eval('v=feval(A_operator,zeros(n,0));','n=-1;')
      if n>-1, J=[J,jj]; end
    end
  else
    if isempty(sigma), sigma=s; end 
    J=[J,jj];  
  end 
end


varg(J)=[];

if n==-1,
  msg1=sprintf('  Cannot find the dimension of ''%s''.  \n',A_operator);
  msg2=sprintf('  Put the dimension n in the parameter list:  \n  like');
  msg3=sprintf('\t\n\n\t jdqr(''%s'',n,..),  \n\n',A_operator);
  msg4=sprintf('  or let\n\n\t n = %s(',A_operator);
  msg5=sprintf('[],''dimension'')\n\n  give n.');
  msg=[msg1,msg2,msg3,msg4,msg5];
  errordlg(msg,'MATRIX')
end

nselect=[]; 
if n<2, return, end

if length(varg) == 1
   nselect=min(n,varg);
elseif length(varg)>1
   if isempty(sigma), sigma=varg(end); varg(end)=[]; end
   nselect=min(n,min(varg));
end

fopts = []; if ~isempty(options), fopts=fieldnames(options); end

if isempty(L_precond) 
  if strmatch('Precond',fopts) 
     L_precond = options.Precond;
  elseif strmatch('L_Precond',fopts)
     L_precond = options.L_Precond;
  end
end

if isempty(U_precond) & strmatch('U_Precond',fopts)
  U_precond = options.U_Precond;
end

if isempty(L_precond), ls_tol  = [0.7,0.49]; end
n=SetPrecond(n); if n<0, return, end

if strmatch('Tol',fopts),   tol = options.Tol;       end

if isempty(nselect), nselect=min(n,5); end

if strmatch('jmin',fopts), jmin=min(n,options.jmin); end 
if strmatch('jmax',fopts), jmax=min(n,options.jmax); end
if jmin<=0, jmin=1; end
if jmax<=jmin, jmax=jmin+1; end

if strmatch('MaxIt',fopts), maxit = abs(options.MaxIt);   end

if strmatch('v0',fopts) ;
 if size(options.v0,2)>0 
  m=size(options.v0,1); 
  V=options.v0(1:min(m,n),1);
  V = [V; ones(n-m,1)];
  nr=norm(V);
  if nr>0, V=V/nr; else,  V = ones(n,1)/sqrt(n) ; end 
 else
  V = ones(n,1)/sqrt(n) ;
 end
else
  V = ones(n,1)/sqrt(n) ;
end


if strmatch('Disp',fopts),     SHOW  = boolean(options.Disp,[SHOW,2]);  end

if SHOW
   ShowChoices(n,nselect,sigma,jmin,jmax,tol,maxit,V,SHOW)
end

return
%-------------------------------------------------------------------
function n=SetPrecond(n)
% finds out how the preconditioners are defined (type_preconditioner)
% and checks consistency of the definitions.
%
% If M is the preconditioner then P*M=L*U. Defaults: L=U=P=I.
%
% type_preconditioner
%       0:   no L
%       1:   L M-file, no U,     L ~= A
%       2:   L M-file, no U,     L == A
%       3:   L M-file, U M-file, L ~= A, U ~= A, L ~=U
%       4:   L M-file, U M-file, L == U
%       5:   L matrix, no U
%       6:   L matrix, U matrix  no P
%       7:   L matrix, U matrix, P matrix

  global A_operator ...
         type_preconditioner ...
         L_precond U_precond P_precond

  % Set type preconditioner
  type_preconditioner=0;
  if isempty(L_precond), return, end

  if ~isempty(U_precond) && ischar(L_precond)~=ischar(U_precond)
    msg=sprintf('  L and U should both be strings or matrices');
    errordlg(msg,'PRECONDITIONER'), n=-1; return
  end
  if ~isempty(P_precond) && (ischar(P_precond) || ischar(L_precond))
    msg=sprintf('  P can be specified only if P, L and U are matrices'); 
    errordlg(msg,'PRECONDITIONER'), n=-1; return
  end  
  tp=1+4*~ischar(L_precond)+2*~isempty(U_precond)+~isempty(P_precond);
  if tp==1, tp = tp + strcmp(L_precond,A_operator); end
  if tp==3, tp = tp + strcmp(L_precond,U_precond); end
  if tp==3 && strcmp(U_precond,A_operator)
    msg1=sprintf('  If L and A use the same M-file,')
    msg2=sprintf('\n  then so should U.'); 
    errordlg([msg1,msg2],'PRECONDITIONER'), n=-1; return
  end
  if tp>5, tp=tp-1; end, type_preconditioner=tp;

  % Check consistency definitions
  if tp<5 && exist(L_precond) ~=2
    msg=sprintf('  Can not find the M-file ''%s.m''  ',L_precond); 
    errordlg(msg,'PRECONDITIONER'), n=-1; return
  end

  ok=1; 
  if tp == 2
    eval('v=feval(A_operator,zeros(n,1),''preconditioner'');','ok=0;')
    if ~ok
       msg='Preconditioner and matrix use the same M-file';
       msg1=sprintf(' %s.   \n',L_precond);
       msg2='Therefore the preconditioner is called';
       msg3=sprintf(' as\n\n\tw=%s(v,''preconditioner'')\n\n',L_precond); 
       msg4='Put this "switch" in the M-file.';
       msg=[msg,msg1,msg2,msg3,msg4];
    end
  end

  if tp == 4 || ~ok
    ok1=1;
    eval('v=feval(L_precond,zeros(n,1),''L'');','ok1=0;')
    eval('v=feval(L_precond,zeros(n,1),''U'');','ok1=0;')
    if ok1 
      type_preconditioner = 4; U_precond = L_precond; ok=1;
    else
      if tp == 4
        msg='L and U use the same M-file';
        msg1=sprintf(' %s.m   \n',L_precond);
        msg2='Therefore L and U are called';
        msg3=sprintf(' as\n\n\tw=%s(v,''L'')',L_precond); 
        msg4=sprintf(' \n\tw=%s(v,''U'')\n\n',L_precond); 
        msg5=sprintf('Check the dimensions and/or\n');
        msg6=sprintf('put this "switch" in %s.m.',L_precond);
        msg=[msg,msg1,msg2,msg3,msg4,msg5,msg6]; 
      end
      errordlg(msg,'PRECONDITIONER'), n=-1; return
    end
  end

  if tp==1 || tp==3
    eval('v=feval(L_precond,zeros(n,1));','ok=0')
    if ~ok
       msg=sprintf('''%s'' should produce %i-vectors',L_precond,n); 
       errordlg(msg,'PRECONDITIONER'), n=-1; return 
    end
  end

  if tp==3
    if exist(U_precond) ~=2
      msg=sprintf('  Can not find the M-file ''%s.m''  ',U_precond);
      errordlg(msg,'PRECONDITIONER'), n=-1; return
    else
      eval('v=feval(U_precond,zeros(n,1));','ok=0')
      if ~ok
        msg=sprintf('''%s'' should produce %i-vectors',U_precond,n);
        errordlg(msg,'PRECONDITIONER'), n=-1; return 
      end
    end
  end
  
  if tp==5
    if min([n,2*n]==size(L_precond)) 
      U_precond=L_precond(:,n+1:2*n); L_precond=L_precond(:,1:n); 
      type_preconditioner=6;
    elseif min([n,3*n]==size(L_precond)) 
      U_precond=L_precond(:,n+1:2*n); P_precond=L_precond(:,2*n+1:3*n);
      L_precond=L_precond(:,1:n); type_preconditioner=7;
    elseif ~min([n,n]==size(L_precond)) 
      msg=sprintf('The preconditioning matrix\n');
      msg2=sprintf('should be %iX%i or %ix%i ([L,U])\n',n,n,n,2*n); 
      msg3=sprintf('or %ix%i ([L,U,P])\n',n,3*n); 
      errordlg([msg,msg2,msg3],'PRECONDITIONER'), n=-1; return
    end
  end

  if tp==6 && ~min([n,n]==size(L_precond) & [n,n]==size(U_precond))
    msg=sprintf('Both L and U should be %iX%i.',n,n); n=-1;
    errordlg(msg,'PRECONDITIONER'), n=-1; return 
  end

  if tp==7 && ~min([n,n]==size(L_precond) & ...
       [n,n]==size(U_precond) & [n,n]==size(P_precond))
    msg=sprintf('L, U, and P should all be %iX%i.',n,n); n=-1;
    errordlg(msg,'PRECONDITIONER'), n=-1; return 
  end

return
%-------------------------------------------------------------------
function x = boolean(x,gamma,string)
%Y = BOOLEAN(X,GAMMA,STRING)
%  GAMMA(1) is the default. 
%  If GAMMA is not specified, GAMMA = 0.
%  STRING is a matrix of accepted strings. 
%  If STRING is not specified STRING = ['no ';'yes']
%  STRING(I,:) and GAMMA(I) are accepted expressions for X 
%  If X=GAMMA(I) then Y=X. If X=STRING(I,:), then Y=GAMMA(I+1).
%  For other values of X, Y=GAMMA(1);

if nargin < 2, gamma=0; end
if nargin < 3, string=strvcat('no','yes'); gamma=[gamma,0,1]; end

if ischar(x)
  i=strmatch(lower(x),string,'exact'); 
  if isempty(i),i=1; else, i=i+1; end, x=gamma(i);
elseif max((gamma-x)==0)
elseif gamma(end) == inf
else, x=gamma(1);
end
  
return
%===========================================================================
%============= OUTPUT FUNCTIONS ============================================
%===========================================================================
function ShowOptions
 fprintf('\n')
 fprintf('PROBLEM\n')
 fprintf('            A: [ square matrix | string ]\n');
 fprintf('      nselect: [ positive integer {5} ]\n\n');
 fprintf('SHIFT\n')
 fprintf('        sigma: [ scalar | row or vector of scalars | {[]} ]\n');

 fprintf('OPTIONS\n');
 fprintf('          Tol: [ positive scalar {1e-8} ]\n');
 fprintf('         Disp: [ yes | {no} ]\n');
 fprintf('         jmin: [ positive integer {nselect+5} ]\n');
 fprintf('         jmax: [ positive integer {jmin+5} ]\n');
 fprintf('        MaxIt: [ positive integer {200} ]\n');
 fprintf('           v0: [ size(A,1) vector of scalars {ones(size(A,1),1)} ]\n');
 fprintf('\nPRECONDITIONER\n');
 fprintf('            M: [ n by n or n by 2n matrix {identity} | string ]\n');
 fprintf('or M=L*U (strings or matrices), or M=P\\(L*U) (matrices) where\n');
 fprintf('            L: [ n by n matrix {identity} | string ]\n');
 fprintf('            U: [ n by n matrix {identity} | string ]\n');
 fprintf('            P: [ n by n {identity} ]\n');
 fprintf('\n')

return
%------------------------------------------------------------------------
function ShowChoices(n,nselect,sigma,...
         jmin,jmax,tol,maxit,V,SHOW)

global A_operator type_preconditioner L_precond U_precond P_precond

  fprintf('\n'),fprintf('PROBLEM\n')
  if ischar(A_operator)
    fprintf('            A: ''%s''\n',A_operator);
  elseif issparse(A_operator)
    fprintf('            A: [%ix%i sparse]\n',n,n);
  else
    fprintf('            A: [%ix%i double]\n',n,n);
  end
    fprintf('    dimension: %i\n',n);
    fprintf('      nselect: %i\n\n',nselect);

  if isreal(sigma)&&isnumeric(sigma)&&~isempty(sigma)
    fprintf('SHIFT\n')
    Str=ShowLambda(sigma);
    fprintf('        sigma: %s',Str)
  end
    fprintf('\n\n')

  fprintf('OPTIONS\n');
    fprintf('          Tol: %g\n',tol);
    fprintf('         Disp: %i\n',SHOW);
    fprintf('         jmin: %i\n',jmin);
    fprintf('         jmax: %i\n',jmax);
    fprintf('        MaxIt: %i\n',maxit);
    fprintf('           v0: [%ix%i double]\n',size(V));

 fprintf('\nPRECONDITIONER\n');
  if type_preconditioner == 0
    fprintf('    No preconditioner\n');
  elseif type_preconditioner == 1
    fprintf('            M: ''%s''\n',L_precond);
  elseif type_preconditioner == 2
    fprintf('            M: ''%s'' [''preconditioner'']\n',L_precond);
  elseif type_preconditioner == 3
    fprintf('            M = L*U where\n');
    fprintf('            L: ''%s''\n',L_precond);
    fprintf('            U: ''%s''\n',U_precond);
  elseif type_preconditioner == 4
    fprintf('    M = L*U where\n');
    fprintf('            L: ''%s'' [''L'']\n',L_precond);
    fprintf('            U: ''%s'' [''U'']\n',U_precond);
  else
    Ls='double'; if issparse(L_precond), Ls='sparse'; end
    if type_preconditioner == 5
      fprintf('%13s: [%ix%i %s]\n','M',n,n,Ls);
    else
      fprintf('            M = L*U where\n');
      fprintf('            L: [%ix%i %s]\n',n,n,Ls);
      Us='double'; if issparse(U_precond), Us='sparse'; end
      fprintf('            U: [%ix%i %s]\n',n,n,Us);
      if type_preconditioner == 7
        fprintf('            P: [%ix%i %s]\n',n,n,Us);
      end
    end
  end

  fprintf('\n\n')

return
%-------------------------------------------------------------------
function varargout=ShowLambda(lambda,kk)

for k=1:size(lambda,1);
  if k>1, Str=[Str,sprintf('\n%15s','')]; else, Str=[]; end
  rlambda=real(lambda(k,1)); ilambda=imag(lambda(k,1));
  Str=[Str,sprintf(' %+11.4e',rlambda)];
  if abs(ilambda)>100*eps*abs(rlambda)
    if ilambda>0
      Str=[Str,sprintf(' + %10.4ei',ilambda)];
    else
      Str=[Str,sprintf(' - %10.4ei',-ilambda)];
    end
  end
end

if nargout == 0
  if nargin == 2
    Str=[sprintf('\nlambda(%i) =',kk),Str];
  else
    Str=[sprintf('\nDetected eigenvalues:\n\n%15s',''),Str];
  end
  fprintf('%s\n',Str)
else
  varargout{1}=Str;
end

return
%===========================================================================
function DispResult(s,nr,gamma)

  if nargin<3, gamma=0; end

  extra='';
  if nr > 100*eps && gamma
     extra='    norm > 100*eps !!! ';
  end

  if gamma<2 || nr>100*eps
     fprintf('\n %35s: %0.5g\t%s',s,nr,extra)
  end

return
%===========================================================================

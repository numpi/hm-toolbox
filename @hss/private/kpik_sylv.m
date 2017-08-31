function [X1,X2,res]=kpik_sylv(A,LA,UA,B,LLB,UB,rhs1,rhs2,m,tol)
%function [X1,X2,res]=kpik_sylv(A,LA,UA,B,LLB,UB,rhs1,rhs2,m,tol)
%
% KPIK algorithm for Sylvester equation  (extended Krylov subspace type method)
%
%   solve    A X + X B' + rhs1 rhs2' = 0    (mind the transposition in B !!)
%                                            rhs1,rhs2 tall matrices
%  assume sep(A,B) >0 
%
%  Author of the matlab code: V. Simoncini
%  Version 1.0
%  Please refer to V.Simoncini for problems using this code.
%  The code uses the function lyap.m of the Control Matlab Toolbox
%
%  Input:
%  A = LA*UA;   LU factorization of A,   size(A) = na x na
%  B = LLB*UB;   LU factorization of B,  size(B) = nb x nb,  nb possibly different from na
%  rhs1, rhs2   factors of right-hand side
%  m = max space dimension, say min(sqrt(size(A)),sqrt(size(B)))
%  tol = max final accuracy (in terms of relative residual)
%  stopping criterion:
%       ||A X1 X2' + X1 X2' B' + rhs1 rhs2'||
%       --------------------------------------  < tol
%                 || rhs1  rhs2'||
%
% Output:
%  X1,X2 = solution factors   X_approx = X1 X2'
%  res  history of residual norms
%
%
%
% If you use this code, please cite the following article:
%
% Tobias Breiten, Valeria Simoncini and Martin Stoll
% Fast iterative solvers for fractional differential equations
% Technical Report, January 2014.
%
%
%THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR 
%IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS 
%FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR 
%COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER 
%IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN 
%CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE. 
%

debug = 0;

nrmb1=norm(rhs1,'fro');  
nrmb2=norm(rhs2,'fro');  
normb=nrmb1*nrmb2;

[n,sh]=size(rhs1);
[nB,sh1]=size(rhs2); if sh~=sh1, fprintf('wrong rhs size\n');return;end

odds=[];
oddsB=[];
res=zeros(m+2,1);

%if norm(A-A',1)<1e-14, kA_max =2 ; else kA_max = m; end
%if norm(B-B',1)<1e-14, kB_max =2 ; else kB_max = m; end

kA_max = 2;
kB_max = 2;

s=2*sh;
rhsA=UA\(LA\rhs1);  
rhsB=UB\(LLB\rhs2);  

% Sequence in A
[U(1:n,1:s),beta]=qr([rhs1,rhsA],0);
ibeta=inv(beta);
beta = beta(1:sh,1:sh); 
H=zeros((m+1)*s,m*s); L=zeros((m+1)*s,m*s);
T=zeros(m*s,m*s);

% Sequence in B'
[W(1:nB,1:s),betaB]=qr([rhs2,rhsB],0);
ibetaB=inv(betaB);
betaB = betaB(1:sh,1:sh); 
HB=zeros((m+1)*s,m*s); LB=zeros((m+1)*s,m*s);
TB=zeros(m*s,m*s);

beta2=beta*betaB';

if debug == 1
    fprintf('      no its      rel res\n')
end

for j=1:m

    jms=(j-1)*s+1;j1s=(j+1)*s;js=j*s;js1=js+1; jsh=(j-1)*s+sh;

% Sequence in A
    Up(1:n,1:sh)   = A*U(1:n,jms:jsh); 
    Up(1:n,sh+1:s) = UA\(LA\U(1:n,jsh+1:js));
%new bases block (modified gram)
    for l=1:2
        k_min=max(1,j-kA_max);
        for kk=k_min:j
            k1=(kk-1)*s+1; k2=kk*s;
            coef= U(1:n,k1:k2)'*Up;
            H(k1:k2,jms:js) = H(k1:k2,jms:js)+ coef; 
            Up = Up - U(:,k1:k2)*coef; 
        end
    end
    if (j<=m)
       [U(1:n,js1:j1s),H(js1:j1s,jms:js)]=qr(Up,0);
       hinv=inv(H(js1:j1s,jms:js));
    end
    I=speye(js+s);
    if (j==1)
      L(1:j*s+sh,(j-1)*sh+1:j*sh) =...
      [ H(1:s+sh,1:sh)/ibeta(1:sh,1:sh), speye(s+sh,sh)/ibeta(1:sh,1:sh)]*ibeta(1:s,sh+1:s);
    else
      L(1:j*s+s,(j-1)*sh+1:j*sh) = L(1:j*s+s,(j-1)*sh+1:j*sh) + H(1:j*s+s,jms:jms-1+sh)*rho;
    end
    odds  = [odds, jms:(jms-1+sh)];   % store the odd block columns
    evens = 1:js; evens(odds)=[];
    T(1:js+s,odds)=H(1:js+s,odds);   %odd columns

    T(1:js+sh,evens)=L(1:js+sh,1:j*sh);   %even columns
    L(1:j*s+s,j*sh+1:(j+1)*sh) = ...
       ( I(1:j*s+s,(js-sh+1):js)- T(1:js+s,1:js)*H(1:js,js-sh+1:js))*hinv(sh+1:s,sh+1:s);
    rho = hinv(1:sh,1:sh)\hinv(1:sh,sh+1:s);

%Sequence in B'
    Wp(1:nB,1:sh)   = B*W(1:nB,jms:jsh); 
    Wp(1:nB,sh+1:s) = UB\(LLB\W(1:nB,jsh+1:js));

%new bases block (modified gram)
    for l=1:2
        k_min=max(1,j-kB_max);
        for kk=k_min:j
            k1=(kk-1)*s+1; k2=kk*s;
            coef= W(1:nB,k1:k2)'*Wp;
            HB(k1:k2,jms:js) = HB(k1:k2,jms:js)+ coef; 
            Wp = Wp - W(:,k1:k2)*coef; 
        end
    end
    if (j<=m)
       [W(1:nB,js1:j1s),HB(js1:j1s,jms:js)]=qr(Wp,0);
       hinvB=inv(HB(js1:j1s,jms:js));
    end
    I=speye(js+s);
    if (j==1),
      LB(1:j*s+sh,(j-1)*sh+1:j*sh) =...
      [ HB(1:s+sh,1:sh)/ibetaB(1:sh,1:sh), speye(s+sh,sh)/ibetaB(1:sh,1:sh)]*ibetaB(1:s,sh+1:s);
    else
      LB(1:j*s+s,(j-1)*sh+1:j*sh) = LB(1:j*s+s,(j-1)*sh+1:j*sh) + HB(1:j*s+s,jms:jms-1+sh)*rhoB;
    end
    oddsB = [oddsB, jms:(jms-1+sh)];   % store the odd block columns
    evens = 1:js; evens(oddsB)=[];
    TB(1:js+s,oddsB)=HB(1:js+s,oddsB);   %odd columns

    TB(1:js+sh,evens)=LB(1:js+sh,1:j*sh);   %even columns
    LB(1:j*s+s,j*sh+1:(j+1)*sh) = ...
       ( I(1:j*s+s,(js-sh+1):js)- TB(1:js+s,1:js)*HB(1:js,js-sh+1:js))*hinvB(sh+1:s,sh+1:s);
    rhoB = hinvB(1:sh,1:sh)\hinvB(1:sh,sh+1:s);

    k=j;

% Solve reduced dim sylvester eqn
    Y = lyap(T(1:js,1:js),TB(1:js,1:js)',eye(k*s,sh)*beta2*eye(k*s,sh)');

    cc  = [H(js1:j1s,js-s+1:js-sh), L(js1:j1s,(j-1)*sh+1:j*sh)];
    ccB = [HB(js1:j1s,js-s+1:js-sh), LB(js1:j1s,(j-1)*sh+1:j*sh)];

    res(k)=sqrt(norm(cc*Y(js-s+1:js,:),'fro')^2+norm(Y(:,js-s+1:js)*ccB','fro')^2)/normb;

    if debug == 1
        disp([k,res(k)])
    end

    if (res(k)<tol), break, end
end

X1=U(1:n,1:js)*Y; X2=W(1:nB,1:js); 
res=res(1:k);

return


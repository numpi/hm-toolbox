function [Uf, Vf] = hm_lanczos(f1,f2,m,n,tol)
%-----------------------------------------------------------------------------------------
% Function that computes a low rank approximation of the matrix A ~ Uf * Vf'
% A is a (m x n)-matrix with some special structure that allows fast matrix-vector product
% f1(x) = Ax
% f2(x) = A'x
%-----------------------------------------------------------------------------------------
if ~exist('tol','var')
	tol = 1e-8;
end
k_max = min([50, m, n]);
U = zeros(m,k_max);
V = zeros(n,k_max);
alpha = zeros(k_max,1); 
betha = zeros(k_max-1,1);
U(:,1) = randn(m,1);
%U(:,1) = f1(eye(n,1));
U(:,1) = U(:,1)/norm(U(:,1));
V(:,1) = f2(U(:,1));
alpha(1) = norm(V(:,1));
V(:,1)= V(:,1)/alpha(1);
omega = 1;
it = 2;
while 1
	omega_old=omega;
	if it > k_max
		U=[U,zeros(m,k_max)];
		V=[V,zeros(n,k_max)];
		alpha=[alpha;zeros(k_max,1)];
		betha=[betha;zeros(k_max-1,1)];
		k_max=2*k_max;
	end
	
	U(:,it) = f1(V(:,it-1)) - alpha(it-1) * U(:,it-1);

	U(:,it)=U(:,it) - U(:,1:it-1)*(U(:,1:it-1)'*U(:,it)); % re-orthogonalization step
	U(:,it)=U(:,it) - U(:,1:it-1)*(U(:,1:it-1)'*U(:,it));

	betha(it-1) = norm(U(:,it));
	U(:,it) = U(:,it)/betha(it-1);
	
	V(:,it) = f2(U(:,it)) - betha(it-1) * V(:,it-1);
	
	V(:,it)=V(:,it) - V(:,1:it-1)*(V(:,1:it-1)'*V(:,it)); % re-orthogonalization step
	V(:,it)=V(:,it) - V(:,1:it-1)*(V(:,1:it-1)'*V(:,it));

	alpha(it) = norm(V(:,it));
	V(:,it) = V(:,it)/alpha(it);
	  [it, betha(it-1),alpha(it), tol]
	if(betha(it-1)<tol && alpha(it)<tol)
		break
	end
	it = it+1;
end
it = it-1;
[alpha(1:it),[0;betha(1:it-1)]];
[Uf,S,Vf]=svd(diag(alpha(1:it))+diag(betha(1:it-1),-1));

Uf=U(:,1:it)*Uf*S;
Vf=V(:,1:it)*Vf;

function x = hss_chol_fact_solve(F,b)
% HSS_CHOL_FACT_SOLV computes the solution of the linear system A X =B
%	       where F contains the generalized Cholesky decomposition of A
%
%
x = b;
for j = 1:length(F.U)
	x(F.ind{j}, :) = F.L{j} \ x(F.ind{j}, :);
	x(F.ind{j}, :) = F.U{j}' * x(F.ind{j}, :);				
end 
x(F.ind{end}, :) = F.L{end}'\(F.L{end} \ x(F.ind{end}, :));
for j = length(F.U):-1:1
		x(F.ind{j}, :) = F.U{j} * x(F.ind{j}, :);
		x(F.ind{j}, :) = F.L{j}' \ x(F.ind{j}, :);				
end 

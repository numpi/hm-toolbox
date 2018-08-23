function p=myadipars_n(a,b,alpha,l0)
%
% function p=myadipars(a,b,alpha,l0)
%
% calculates the optimal ADI shiftparameters (for equations where
% Matrix A is stable and symmetric) according to Jing-Rebecca Li
% and Jakob Whites "Low Rank Solution of Lyapunov equation" (which
% gives an overview of Wachspress`s method form e.g. "The ADI model Problem"
%
% p   is the array of shift parameters
%
% a   is assumed to be the absolute value of the smallest magnitude
%     eigenvalue of the Systemmatrix A
%
% b   is assumed to be the absolute value of the largest magnitude eigenvalue
%
% alpha is the arctan of the maximum of abs(imag(lamba))/abs(real(lambda))
%       for all lambda in the spectrum of A
%
% l0  is the number of desired shift parameters that should be calculated
%
if (alpha==0)
    kprime=a/b;
else
    c2 = 2/(1+(a/b+b/a)/2);
    m = 2*cos(alpha)*cos(alpha)/c2 -1;
    if (m<1)
        error(['Shift parameters would be complex, function not aplicable, ' ...
            'aborting!']);
        
        %
        % FIX ME: if m<1 parameter become complex! switch back to the
        % heuristics by Thilo or complex parameters suggested by
        % Wachspress.
        %
        % ALA2006 -> also test switching to Krylov projection based
        % method see V.Simoncini
        %
    end
    kprime = 1/(m+sqrt(m^2-1));
end
k=min(1-eps,sqrt(1-kprime^2));%this is a workaround for the case
%k=1 that works for Peters Model
%reduction problems not really
%nice but it works great.

%TODO: check the computation of k, kprime to avoid roundoff errors
%and probably replace the hack above.

[K,E]=ellip(k,pi/2);
if (alpha==0)
    [v,E]=ellip(kprime,pi/2);
else
    [v,E]=ellip(kprime,asin(sqrt(a/(b*kprime))));
end
J=l0;

p=ones(J,1);
for i=1:J
    p(i)=-sqrt(a*b/kprime)*dn((i-0.5)*K/J,k);
    %here we have the choice to take the
    %matlab function ellipj or my own
    %one dn. the later can be proted to
    %FORTRAN or C Code very easily
end
p=-p;

function p=myadipars(a,b,alpha,TOL)
%
% function p=myadipars(a,b,alpha,TOL)
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
% TOL is the epsilon1 of the above paper. The smaller the better
%     the approximation but the more parameters are calculated
%

%%% Davide's comment
%%% also the notation follows the paper by Jing-Rebecca Li
%%% and Jakob Whites "Low Rank Solution of Lyapunov equation", pp. 264-265

if (alpha==0)
    kprime=a/b;
else
    cos2_beta = 2/(1+(a/b+b/a)/2);
    m = 2*cos(alpha)*cos(alpha)/cos2_beta -1;
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
if (alpha==0) %%% Davide's comment: just because if alpha==0 then kprime=a/b and therefore asin(sqrt(a/(b*kprime)))=pi/2
    [v,E]=ellip(kprime,pi/2);
else
    [v,E]=ellip(kprime,asin(sqrt(a/(b*kprime))));
end
%%% Davide's comment
%%% REMARK: given a tolerance TOL we can known in advanced the number of shifts needed to get that accuracy!
J=ceil(K/(2*v*pi)*log(4/TOL));

p=ones(J,1);
for i=1:J
    p(i)=-sqrt(a*b/kprime)*dn((i-0.5)*K/J,k);
    %here we have the choice to take the
    %matlab function ellipj or my own
    %one dn. the later can be proted to
    %FORTRAN or C Code very easily
end

%%%% we want p to be positive
p=-p;

function out=dn(u,k)
%function out=dn(u,k) calculates the value of the elliptic
%function dn (see Abramowitz/Stegun Handbook of mathematical
%functions '65)
%

a(1)=1;
c(1)=k;
b(1)=min(1-eps,sqrt(1-k*k));
i=1;
while abs(c(i))>eps
    i=i+1;
    a(i)=(a(i-1)+b(i-1))/2;
    b(i)=sqrt(a(i-1)*b(i-1));
    c(i)=(a(i-1)-b(i-1))/2;
end
phi1=(2.^(i-1)).*a(i).*u;%here 2^(i-1) and not 2^i as in
%Abramowitz/Stegun because counting
%starts at 1 not at 0 like in the book
phi0=0;
for j=i:-1:2
    if (j<i)
        phi1=phi0;
    end;
    phi0=(phi1+asin(c(j)*sin(rem(phi1,2*pi))/a(j)))/2;
end
arg=1-k*k*sin(rem(phi0,2*pi))^2;
if (arg<.1)
    out=sqrt(arg);
else
    out=cos(rem(phi0,2*pi))/cos(phi1-phi0);
end

%the last two are both representations found in the
%Abramowitz/Stegun book. if arg is close to zero the cosine
%version should be better to avoid numerical inexactness resulting
%from the substraction.
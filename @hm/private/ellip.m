function [F,E]=ellip(hk,phi);
%  function [F,E]=ellip(hk,phi);
%  Computes complete and incomplete elliptic integrals F(k,phi) and E(k,phi)
%       Input  : hk  --- Modulus k ( 0 < k < 1 )
%                phi --- Argument 
%       Output : F   --- F(k,phi)
%                E   --- E(k,phi)
%       ==================================================

g=0.0;
a0=1.0;
b0=min(1-eps,sqrt(1.0-hk.*hk));
d0=phi;
r=hk.*hk;
if (hk == 1.0&phi == pi/2) ;
  F=1.0e+300;
  E=1.0;
elseif (hk == 1.0);
  F=log((1.0+sin(d0))./cos(d0));
  E=sin(d0);
else;
  fac=1.0;
  for  n=1:40;
    a=(a0+b0)./2.0;
    b=sqrt(a0.*b0);
    c=(a0-b0)./2.0;
    fac=2.0.*fac;
    r=r+fac.*c.*c;
    if (phi ~= pi/2) ;
      d=d0+atan((b0./a0).*tan(d0));
      g=g+c.*sin(d);
      d0=d+pi.*fix(d./pi+.5);
    end;
    a0=a;
    b0=b;
    if (c < 1.0e-15) break; end;
  end;
  ck=pi./(2.0.*a);
  ce=pi.*(2.0-r)./(4.0.*a);
  if (phi == pi/2) ;
    F=ck;
    E=ce;
  else
    F=d0./(fac.*a);
    E=F.*ce./ck+g;
  end;
end;
return;

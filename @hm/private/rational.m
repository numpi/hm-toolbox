function [R,P,K]=rational(nu,type)
%function [R,P,K]=rational(nu,type)
%
%      B(s)       R(1)       R(2)             R(n)
%       ----  =  -------- + -------- + ... + -------- + K(s)
%      A(s)     s - P(1)   s - P(2)         s - P(n)
%
% type =1 cheb, else  Pade'
if (type==1)  %Cheb  (only nu=10 at the moment)
    [B,A]=cheby(nu);
else           %Pade
    [B,A]=coeff(nu);
end

[R,P,K]=residue(B(nu+1:-1:1),A(nu+1:-1:1));
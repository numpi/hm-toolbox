function [X,Y] = sda(E,F,G,P)
% [X,Y]=SDA(E,F,G,P) applies the Structured Doubling Algorithm starting 
% from E, F, G, P
% E, F, G, P: initial matrices
% X: limit of the sequence (P_k)
% Y: limit of the sequence (G_k)
tol = 1e-13;
kmax = 30;
err = 1;
k = 0;
n = size(G,1); m = size(P,1);
while err > tol && k < kmax
    IGP = eye(n) - G*P;
    IPG = eye(m) - P*G;
    E1 = E/IGP;
    F1 = F/IPG;
    G = G + E1*G*F;
    P = P + F1*P*E;
    E = E1*E;
    F = F1*F;
    err = min(norm(E,1),norm(F,1));
    k = k + 1;
end
success = 1;
X = P; Y = G;
if k == kmax 
    disp('SDA:: Warning: reached the maximum number of iterations')
end

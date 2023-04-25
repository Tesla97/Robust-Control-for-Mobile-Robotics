function [K] = controlloreH2(A,B1,B2,C2,D22,n,m,abilitazione,a,theta,r)
% Definizione matrici incognite
X     = sdpvar(n,n);
W     = sdpvar(m,n);
Z     = sdpvar(1,1);
gamma = sdpvar(1,1);
% LMIs H2 minimizzazione
F1 = ([(A*X+B1*W)+(A*X+B1*W)'  B2;
                 B2'          -eye(1)]<= 0);
F2 = ([      Z  (C2*X+D22*W);
       (C2*X+D22*W)'  X    ]>= 0);
F3 = (trace(Z) <= gamma);
F4 = ([gamma] >= 0.02);
% D-stabilità LMIs
F5 = ([(A*X+B1*W)+(A*X+B1*W)'+2*a*X]<=0);
F6 = ([-r*X  (A*X+B1*W);(A*X+B1*W)'  -r*X]<=0);
F7 = ([(A*X+B1*W)*sin(theta)+(A*X+B1*W)'*sin(theta)  (A*X+B1*W)*cos(theta)-(A*X+B1*W)'*cos(theta);
       (A*X+B1*W)'*cos(theta)-(A*X+B1*W)*cos(theta)  (A*X+B1*W)*sin(theta)+(A*X+B1*W)'*sin(theta)]<= 0);
% scelta opzione D-stabilità
if abilitazione
    F = F1+F2+F3+F4+F5+F6+F7;
else
    F = F1+F2+F3+F4;
end
% ottimizzazione
ops= sdpsettings('solver','sedumi');
optimize(F,[gamma],ops);
% guadagno
X = double(X);
W = double(W);
K = W / X;
end
function [K] = controlloreHinfinito(A,B1,B2,C1,D1,D2,n,m,abilitazione,a,theta,r)
% Definizione Matrici incognite
X     = sdpvar(n,n);
W     = sdpvar(m,n);
gamma = sdpvar(1,1);
% LMIs Hinf minimizzazione
F1    = ([X] >= 0);
F2    = ([(A*X+B1*W)+(A*X+B1*W)'    B2        (C1*X+D1*W)'   ;
                    B2'          -gamma*eye(1)        D2'    ;
                (C1*X+D1*W)          D2          -gamma*eye(1)]<= 0);
F3    = ([gamma] >= 0.02);
% D-stabilità LMIs
F4 = ([(A*X+B1*W)+(A*X+B1*W)'+2*a*X]<=0);
F5 = ([-r*X  (A*X+B1*W);(A*X+B1*W)'  -r*X]<=0);
F6 = ([(A*X+B1*W)*sin(theta)+(A*X+B1*W)'*sin(theta)  (A*X+B1*W)*cos(theta)-(A*X+B1*W)'*cos(theta);
       (A*X+B1*W)'*cos(theta)-(A*X+B1*W)*cos(theta)  (A*X+B1*W)*sin(theta)+(A*X+B1*W)'*sin(theta)]<= 0);
% opzione D-stabilità
if abilitazione
    F = F1+F2+F3+F4+F5+F6;
else
    F = F1+F2+F3;
end
%ottimizzazione
ops   = sdpsettings('solver','sedumi');
optimize(F,[gamma],ops);
% guadagano 
X      = double(X);
W      = double(W);
K      = W / X;
end
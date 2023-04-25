function [K] = controlloreL1(A,B1,B2,C3,D31,D32,n,m,abilitazione,a,theta,r)
% lambda 
lambda = 1;
% Definizione Variabili
X     = sdpvar(n,n);
W     = sdpvar(m,n);
eta   = sdpvar(1,1);
gamma = sdpvar(1,1);
% LMIs L1 minimizzazione
F1    = ([(A*X+B1*W)+(A*X+B1*W)'+lambda*X        B2;
                         B2'                -eta*eye(1)]<=0);             
F2    = ([ lambda*X        zeros(4,1)         (C3*X+D32*W)';
          zeros(1,4)   (gamma-eta)*eye(1)          D31';
        (C3*X+D32*W)         D31                gamma*eye(1)]>= 0);   
F3    = ([X] >= 0);
F4    = ([gamma] >= 0.02);
% D-stabilità LMIs
F5    = ([eta] >= 0);
F6 = ([(A*X+B1*W)+(A*X+B1*W)'+2*a*X]<=0);
F7 = ([-r*X  (A*X+B1*W);(A*X+B1*W)'  -r*X]<=0);
F8 = ([(A*X+B1*W)*sin(theta)+(A*X+B1*W)'*sin(theta)  (A*X+B1*W)*cos(theta)-(A*X+B1*W)'*cos(theta);
       (A*X+B1*W)'*cos(theta)-(A*X+B1*W)*cos(theta)  (A*X+B1*W)*sin(theta)+(A*X+B1*W)'*sin(theta)]<= 0);
% scelta opzione D-stabilità
if abilitazione
    F = F1+F2+F3+F4+F5+F6+F7+F8;
else
    F = F1+F2+F3+F4+F5;
end
% ottimizzazione
ops   = sdpsettings('solver','sedumi');
optimize(F,[gamma],ops);
% Guadagno
X     = double(X);
W     = double(W);
K     = W / X;
end
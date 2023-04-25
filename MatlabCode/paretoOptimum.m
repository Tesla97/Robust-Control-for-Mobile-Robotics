function [K] = paretoOptimum(A,B1,B2,C,D,C1,C2,n,m)
% Definizione Matrici Incognite
X        = sdpvar(n,n);
W        = sdpvar(m,n);
Z        = sdpvar(1,1);
gammaInf = sdpvar(1,1);
gamma2   = sdpvar(1,1);
% Definizione Coeffienti alpha e beta
alpha    = 50;
beta     = 25;
% Definizione LMIs Minimizzazione Norma H2
D22      = 0;
F1       = ([(A*X+B1*W)+(A*X+B1*W)'  B2;
                      B2'          -eye(1)]<= 0);
F2       = ([      Z  (C2*X+D22*W);
             (C2*X+D22*W)'  X    ]>= 0);
F3       = (trace(Z) <= gamma2);
% Definizione LMIs Minimizzazione Norm Hinf
D1       = 0;
D2       = 0;
F4       = ([(A*X+B1*W)+(A*X+B1*W)'    B2        (C1*X+D1*W)'   ;
                    B2'          -gammaInf*eye(1)        D2'    ;
                (C1*X+D1*W)          D2          -gammaInf*eye(1)]<= 0);
F5       = ([gamma2 >= 0.02]);
F6       = ([gammaInf >= 0.02]);
% Risoluzione Problema
F        = F1+F2+F3+F4+F5+F6;
ops      = sdpsettings('solver','sedumi');
optimize(F,[alpha*gammaInf + beta*gamma2],ops);
% guadagano 
X        = double(X);
W        = double(W);
K        = W / X;
end
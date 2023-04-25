function [K] = controlloreRobustoH2(A_f,B1_f,B2_f,lmin,lmax,Ipmin,Ipmax,C2,D22,n,m)
% Definizione Vertici Politopo

%Vertice 1
A1 = A_f (lmin,Ipmin);
B11= B1_f(lmin,Ipmin);
B21= B2_f(lmin,Ipmin);
%Vertice 2
A2 = A_f (lmax,Ipmin);
B12= B1_f(lmax,Ipmin);
B22= B2_f(lmax,Ipmin);
%Vertice 3
A3 = A_f (lmin,Ipmax);
B13= B1_f(lmin,Ipmax);
B23= B2_f(lmin,Ipmax);
%Vertice 4
A4 = A_f (lmax,Ipmax);
B14= B1_f(lmax,Ipmax);
B24= B2_f(lmax,Ipmax);

% Definizione matrici incognite
X     = sdpvar(n,n);
W     = sdpvar(m,n);
Z     = sdpvar(1,1);
gamma = sdpvar(1,1);

% LMIs
F1    = ([(A1*X+B11*W)+(A1*X+B11*W)'  B21;
                 B21'           -eye(1)]<= 0);
F2    = ([(A2*X+B12*W)+(A2*X+B12*W)'  B22;
                 B22'           -eye(1)]<= 0);
F3    = ([(A3*X+B13*W)+(A3*X+B13*W)'  B23;
                 B23'           -eye(1)]<= 0);
F4    = ([(A4*X+B14*W)+(A4*X+B14*W)'  B24;
                 B24'           -eye(1)]<= 0);
F5    = ([      Z  (C2*X+D22*W);
           (C2*X+D22*W)'  X    ]>= 0);
F6    = (trace(Z) <= gamma);
F7    = ([gamma] >= 0.02);

% Ottimizzazione
F     = F1+F2+F3+F4+F5+F6+F7;
ops   = sdpsettings('solver','sedumi');
optimize(F,[gamma],ops);

% Guadagno
X     = double(X);
W     = double(W);
K     = W / X;
end
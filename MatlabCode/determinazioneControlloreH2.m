function [Ak,Bk,Ck,Dk] = determinazioneControlloreH2(A,B1,B2,C,D,C2,D22,n,m,p)
% Definizione Matrici Incognite
X     = sdpvar(n,n);
Y     = sdpvar(n,n);
Akc   = sdpvar(n,n);
Bkc   = sdpvar(n,p);
Ckc   = sdpvar(m,n);
Dkc   = sdpvar(m,p);
Q     = sdpvar(1,1);

gamma = 0.02;

%LMIs
F1    = ([(A*Y+B1*Ckc)+(A*Y+B1*Ckc)'    Akc'+(A+B1*Dkc*C)       (B2+B1*Dkc*D);
           Akc+(A+B1*Dkc*C)'         (X*A+Bkc*C)+(X*A+Bkc*C)'   (X*B2+Bkc*D) ;
           (B2+B1*Dkc*D)'               (X*B2+Bkc*D)'             -eye(1)]<=0);
F2    = ([       Q         (C2*Y+D22*Ckc)  (C2+D22*Dkc*C);
          (C2*Y+D22*Ckc)'        Y            eye(4)     ;
          (C2+D22*Dkc*C)'       eye(4)          X       ]>=0);
F3    = ([trace(Q) <= gamma]);
F4    = ([Y eye(4);eye(4) X]>= 0);
F     = F1+F2+F3+F4;
%ottimizzazione
ops = sdpsettings('solver','sedumi');
solvesdp(F,gamma,ops);
% Fattorizzazione LU
X     = double(X);
Y     = double(Y);
P     = eye(4) - X*Y;
[L,U] = lu(P);

% Definizione Matrici Di Inversione
M     = U';
N     = L ;

% Definizione Controllore Dinamico
Akc   = double(Akc);
Bkc   = double(Bkc);
Ckc   = double(Ckc);
Dkc   = double(Dkc);

% Inversione
Dk    = Dkc;
Ck    = (Ckc - Dk*C*Y)*inv(M');
Bk    = inv(N)*(Bkc - X*B1*Dk);
Ak    = inv(N)*(Akc-N*Bk*C*Y-X*B1*Ck*M'-X*(A+B1*Dk*C)*Y)*inv(M');
end
           
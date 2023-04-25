function [Ak,Bk,Ck,Dk] = determinaControlloreHinf(A,B1,B2,C,D,C1,D12,D11,n,m,p)
% Definizione Matrici Incognite
X     = sdpvar(n,n);
Y     = sdpvar(n,n);
Akc   = sdpvar(n,n);
Bkc   = sdpvar(n,p);
Ckc   = sdpvar(m,n);
Dkc   = sdpvar(m,p);

gamma = sdpvar(1,1);

% LMIs
F1  = ([(A*Y+B1*Ckc)+(A*Y+B1*Ckc)'      Akc'+(A+B1*Dkc*C)     (B2+B1*Dkc*D)   (C1*Y+D12*Ckc)' ;
         Akc+(A+B1*Dkc*C)'          (X*A+Bkc*C)+(X*A+Bkc*C)'   (X*B2+Bkc*D)   (C1+D12*Dkc*C)' ;
         (B2+B1*Dkc*D)'                  (X*B2+Bkc*D)'         -gamma*eye(1)  (D11+D12*Dkc*D)';
         (C1*Y+D12*Ckc)                 (C1+D12*Dkc*C)        (D11+D12*Dkc*D)   -gamma*eye(1)]<= 0);
F2  = ([Y eye(4);eye(4) X]>=0);
F3  = ([gamma] >= 0.02);
F   = F1+F2+F3;
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
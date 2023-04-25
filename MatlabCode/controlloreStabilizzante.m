function [K] = controlloreStabilizzante(A,B,n,m)
% Definizione Matrici Incognite
X   = sdpvar(n,n);
W   = sdpvar(m,n);
% Definizione LMIs
F1  =([X] >= 0);
F2  =([(A*X+B*W)+(A*X+B*W)']<= 0);
% Risoluzione LMIs
F   = F1+F2;
ops = sdpsettings('solver','sedumi');
solvesdp(F,[],ops);
% Guadagno K
X   = double(X);
W   = double(W);
K   = W / X;
end
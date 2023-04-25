function [K] = controlloreRobustoStabilizzante(A_f,B1_f,lmin,lmax,Ipmin,Ipmax,n,m)
% Definizione Vertici Politopo

%Vertice 1
A1 = A_f(lmin,Ipmin);
B11= B1_f(lmin,Ipmin);
%Vertice 2
A2 = A_f(lmax,Ipmin);
B12= B1_f(lmax,Ipmin);
%Vertice 3
A3 = A_f(lmin,Ipmax);
B13= B1_f(lmin,Ipmax);
%Vertice 4
A4 = A_f(lmax,Ipmax);
B14= B1_f(lmax,Ipmax);

% Definizione Matrici Incognite
X  = sdpvar(n,n);
W  = sdpvar(m,n);

% LMIs
F1 = ([(A1*X+B11*W)+(A1*X+B11*W)']<= 0);
F2 = ([(A2*X+B12*W)+(A2*X+B12*W)']<= 0);
F3 = ([(A3*X+B13*W)+(A3*X+B13*W)']<= 0);
F4 = ([(A4*X+B14*W)+(A4*X+B14*W)']<= 0);
F5 = ([X] >= 0);

% Risoluzione
F  = F1+F2+F3+F4+F5;
ops= sdpsettings('solver','sedumi');
solvesdp(F);

% Guadagno
X  = double(X);
W  = double(W);
K  = W / X;
end
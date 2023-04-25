function [K] = controlloreDstabilita(A,B,n,m,a,theta,r)
% Definizione Matrici Incognite
X  = sdpvar(n,n);
W  = sdpvar(m,n);
% Definizione LMIs
F1 = ([(A*X+B*W)+(A*X+B*W)'+2*a*X]<=0);
F2 = ([-r*X  (A*X+B*W);(A*X+B*W)'  -r*X]<=0);
F3 = ([(A*X+B*W)*sin(theta)+(A*X+B*W)'*sin(theta)  (A*X+B*W)*cos(theta)-(A*X+B*W)'*cos(theta);
       (A*X+B*W)'*cos(theta)-(A*X+B*W)*cos(theta)  (A*X+B*W)*sin(theta)+(A*X+B*W)'*sin(theta)]<= 0);
% Risoluzione LMIs
F  = F1+F2+F3;
ops = sdpsettings('solver','sedumi');
solvesdp(F,[],ops);
% guadagno
X  = double(X);
W  = double(W);
K  = W / X;
end
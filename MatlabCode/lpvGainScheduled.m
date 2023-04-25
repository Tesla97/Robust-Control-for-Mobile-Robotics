%% Sistemi LPV e Controllo Gain-Scheduled
clear all
close all
clc
%--------------------------------------------------------------------------
% Definizione Parametri Modello
M  = 3;                                       % Massa Cart 
m  = 10.8;                                    % Massa Pendolo al CG 
l  = 0.22;                                    % Altezza CG dal Cart 
g  = 9.81;                                    % Gravità
Wp = 0.20;                                    % Larghezza Piatto Pendolo
Hp = 0.03;                                    % Altezza Piatto Pendolo
Ip = (1/12)*(Wp^2 + Hp^2)*m;                  % Inerzia Rispetto a Yaw
mu = 0.90;                                    % Attrito Pneumatico - Strada

% Definizione Stato Iniziale Sistema
x0 = [0 0 0 0]';

% Parametri Simulazione
ts = 15;                                       % Tempo Di Simulazione
Ts = 0.01;                                     % Tempo Di Campionamento
%--------------------------------------------------------------------------
% Definizione Modello In Stato Del Sistema

% Sistema Lineare Tempo Invariante
% xd = A x + B1 u + B2 w
%  y = C x + D u

% Utilizzo Funzioni Inline Per Il Caso "Modello Incerto"

% Fattori In Comune A Denominatore Nel Modello
mcm1  =            (M+m);
mcm2  = @(l,Ip)    ((m*l^2 + Ip));
com   = @(l,Ip)    ((m*m*l*l)/(mcm1*mcm2(l,Ip)));
mcm3  = @(l,Ip)    (1 - com(l,Ip));

% Definizione Elementi Matrici (A,B,C,D) Del Sistema
A_2_1 = @(l,Ip)    ((m*g*l)/(mcm2(l,Ip)*mcm3(l,Ip)));
A_2_4 = @(l,Ip)    (-1*(m*l*mu)/(mcm1*mcm2(l,Ip)*mcm3(l,Ip)));
A_4_1 = @(l,Ip)    ((m*m*l*l*g)/(mcm1*mcm2(l,Ip)*mcm3(l,Ip)));
A_4_4 = @(l,Ip)    (-mu/(mcm1*mcm3(l,Ip)));

B_1_2 = @(l,Ip)    ((m*l)/(mcm1*mcm2(l,Ip)*mcm3(l,Ip)));
B_1_4 = @(l,Ip)    (1/(mcm1*mcm3(l,Ip)));

B_2_2 = @(l,Ip)    ((l/(mcm2(l,Ip)*mcm3(l,Ip))) - ((m*l)/(mcm1*mcm2(l,Ip)*mcm3(l,Ip))));
B_2_4 = @(l,Ip)    (((m*l*l)/(mcm2(l,Ip)*mcm1*mcm3(l,Ip))) - (1/(mcm1*mcm3(l,Ip))));

% Definizione Matrici Modello Linearizzato
A_f  = @(l,Ip) [    0             1         0            0      ;
               A_2_1(l,Ip)      0         0      A_2_4(l,Ip)  ;
                  0             0         0            1      ;
               A_4_1(l,Ip)      0         0       A_4_4(l,Ip)];
    
B1_f = @(l,Ip) [ 0  B_1_2(l,Ip) 0 B_1_4(l,Ip)]';

B2_f = @(l,Ip) [ 0  B_2_2(l,Ip) 0 B_2_4(l,Ip)]';

A  = A_f(l,Ip);
B1 = B1_f(l,Ip);
B2 = B2_f(l,Ip);

C  = [ 1   0   0   0; 0   0   1   0];

D  = [0 0]';

d  = 25;
%--------------------------------------------------------------------------
% Definizione Range Variabilità Parametri Incerti
lmin = 0.005;
lmax = 0.35;
Ipmin= 0.005;
Ipmax= 0.35;
%--------------------------------------------------------------------------
% Vertici Politopo -> Incertezza Parametrica Politopica

% Vertice 1
AV1  = A_f(lmax,Ipmax);
BV1  = B1_f(lmax,Ipmax);
% Vertice 2
AV2  = A_f(lmax,Ipmin);
BV2  = B1_f(lmax,Ipmin);
% Vertice 3
AV3  = A_f(lmin,Ipmax);
BV3  = B1_f(lmin,Ipmax);
% Vertice 4
AV4  = A_f(lmin,Ipmin);
BV4  = B1_f(lmin,Ipmin);
%--------------------------------------------------------------------------
% Descrizione Politopica 
a    = (l-lmin)/(lmax-lmin);
b    = (Ip-Ipmin)/(Ipmax-Ipmin);
p1   = a*b;
p2   = a*(1-b);
p3   = (1-a)*b;
p4   = (1-a)*(1-b);
%--------------------------------------------------------------------------
% LMIs

% Matrici Incognite
X  = sdpvar(4,4);
W1 = sdpvar(1,4);
W2 = sdpvar(1,4);
W3 = sdpvar(1,4);
W4 = sdpvar(1,4);

% LMIs
F1  = ([(AV1*X+BV1*W1)+(AV1*X+BV1*W1)']<=0);
F2  = ([(AV1*X+BV1*W2)+(AV1*X+BV1*W2)']<=0);
F3  = ([(AV1*X+BV1*W3)+(AV1*X+BV1*W3)']<=0);
F4  = ([(AV1*X+BV1*W4)+(AV1*X+BV1*W4)']<=0);

F5  = ([(AV2*X+BV2*W1)+(AV2*X+BV2*W1)']<=0);
F6  = ([(AV2*X+BV2*W2)+(AV2*X+BV2*W2)']<=0);
F7  = ([(AV2*X+BV2*W3)+(AV2*X+BV2*W3)']<=0);
F8  = ([(AV2*X+BV2*W4)+(AV2*X+BV2*W4)']<=0);

F9  = ([(AV3*X+BV3*W1)+(AV3*X+BV3*W1)']<=0);
F10 = ([(AV3*X+BV3*W2)+(AV3*X+BV3*W2)']<=0);
F11 = ([(AV3*X+BV3*W3)+(AV3*X+BV3*W3)']<=0);
F12 = ([(AV3*X+BV3*W4)+(AV3*X+BV3*W4)']<=0);

F13 = ([(AV4*X+BV4*W1)+(AV4*X+BV4*W1)']<=0);
F14 = ([(AV4*X+BV4*W2)+(AV4*X+BV4*W2)']<=0);
F15 = ([(AV4*X+BV4*W3)+(AV4*X+BV4*W3)']<=0);
F16 = ([(AV4*X+BV4*W4)+(AV4*X+BV4*W4)']<=0);

F17 = ([X]>= 0);

% risoluzione LMIs
F   = F1+F2+F3+F4+F5+F6+F7+F8+F9+F10+F11+F12+F13+F14+F15+F16+F17;
ops = sdpsettings('solver','sedumi');
solvesdp(F,[],ops);

% guadagni
X   = double(X);
W1  = double(W1);
W2  = double(W2);
W3  = double(W3);
W4  = double(W4);
K1  = W1 / X;
K2  = W2 / X;
K3  = W3 / X;
K4  = W4 / X;

%--------------------------------------------------------------------------
% Progettazione stimatore stato (Luemberger)
p = [-10 -11 -12 -13];
H = place(A',C',p)';
As= A - H*C;
Bs=[B1 H];
Cs=[1 0 0 0;0 1 0 0;0 0 1 0;0 0 0 1];
Ds=[0 0 0;0 0 0;0 0 0;0 0 0];

% Gain Scheduled Controller
KLPV   =  p1*K1+p2*K2+p3*K3+p4*K4;





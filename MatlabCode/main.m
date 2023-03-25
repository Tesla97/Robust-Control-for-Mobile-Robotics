%% Codice Matlab Modulo 2 - Nicola Corea 235279 -  Lane Keeping System
clear all
close all
clc
%--------------------------------------------------------------------------
% Definizione Parametri Del Modello
m     = 1573;
Iz    = 2873;
lf    = 1.1;
lr    = 1.58;
Caf   = 80000;
Car   = 80000;
Vx    = 30;
R     = 1000;
L     = lf+lr;
Psidd = (Vx/R);
Ts    = 0.01;

% Definizione Modello In Stato
A_2_2 = (-1)*((2*Caf+2*Car)/(m*Vx));
A_2_3 = ((2*Caf+2*Car)/m);
A_2_4 = ((-2*Caf*lf + 2*Car*lr)/(m*Vx));
A_4_2 = (-1)*((2*Caf*lf-2*Car*lr)/(Iz*Vx));
A_4_3 = ((2*Caf*lf-2*Car*lr)/Iz);
A_4_4 = (-1)*((2*Caf*lf^2 + 2*Car*lr^2)/(Iz*Vx));

B_1_2 = ((2*Caf)/m);
B_1_4 = ((2*Caf*lf)/Iz);

B_2_2 = (-1)*((2*Caf*lf-2*Car*lr)/(m*Vx)) - Vx;
B_2_4 = (-1)*((2*Caf*lf^2 + 2*Car*lr^2)/(Iz*Vx));

% Definizione Matrici
A  = [  0   ,   1   ,   0   ,   0   ;
        0   , A_2_2 , A_2_3 , A_2_4 ;
        0   ,   0   ,   0   ,   1   ;
        0   , A_4_2 , A_4_3 , A_4_4 ];
    
B1 = [0;B_1_2;0;B_1_4];
B2 = [0;B_2_2;0;B_2_4];
C  = [1 0 0 0;0 1 0 0;0 0 1 0;0 0 0 1];
D  = [0 0;0 0;0 0;0 0];

% Stabilizzazione Mediante Pole Placement
% p  = [-5 -6 -7 -8];
% K  = place(A,B1,p);

%Definizione Obiettivo
C1 = [1 0 0 0;0 0 1 0];
D1 = [0 0]';
D2 = [0 1]';

% Problema ottimizzazione
X      = sdpvar(4,4);
W      = sdpvar(1,4);
gamma  = sdpvar(1,1);

% LMIs
F1 = ([X] >= 0);
F2 = ([gamma] >= 0.00001);
F3 = ([(A*X+B1*W)+(A*X+B1*W)'      B2           (C1*X+D1*W)';
                 B2'          -gamma*eye(1)          D2'    ;
              (C1*X+D1*W)          D2           -gamma*eye(2)]<=0);

% Ottimizzazione
F  = F1+F2+F3;
ops= sdpsettings('solver','sedumi');
optimize(F,[gamma],ops);

% Guadagno
X  = double(X);
W  = double(W);
K  = W / X;
K  = -K;
% 
% 
% % Steady State Component
Kv = ((lr*m)/(2*Caf*L))-((lf*m)/(2*Car*L));
% ay = (Vx^2)/R;
% e2r= (-1)*(lr/R) + (lf/(2*Car*L))*((m*Vx*Vx)/R);
% df = (L/R)+Kv*ay+K(3)*e2r;

% Generate Reference
sim("dataSchemeSimulink.slx");
psiDes = ans.psiDesDot;
tsim   = ans.t;
psiDesDot = [tsim psiDes];

figure(1);
plot(tsim,psiDes);
grid;
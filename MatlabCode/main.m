%% Progetto Controllo Dei Veicoli - Nicola Corea - 235279
clear all
close all
clc
%--------------------------------------------------------------------------
% Selezionare Tipo Di Controllore
% 1) controllerType = 1  -> controllore stabilizzante
% 2) controllerType = 2  -> controllore D-stabilità
% 3) controllerType = 3  -> controller minimizzazione norma Hinf  Statico
% 4) controllerType = 4  -> controller minimizzazione norma H2    Statico
% 5) controllerType = 5  -> controller minimizzazione norma L1    Statico
% 6) controllerType = 6  -> confronto prestazioni Hinf , H2 , L1  Statici
% 7) controllerType = 7  -> controller minimizzazione norma Hinf  Dinamico
% 8) controllerType = 8  -> controller minimizzazione norma H2    Dinamico
% 9) controllerType = 9  -> confronto prestazioni Hinf , H2       Dinamici
% 10)controllerType = 10 -> controllore stabilizzante robusto
% 11)controllerType = 11 -> controllore robusto Hinf
% 12)controllerType = 12 -> controllore robusto H2
% 13)controllerType = 13 -> confronto controllori robusti
% 14)controllerType = 14 -> Pareto Optimum Hinf - H2 multiobjective problem
% otherwise              -> pole placement

% Seleziona Tipo Controllore
controllerType = 1;           

% Seleziona D-stabilità (1-> D-stabilità attiva  | 0 -> No)
abilitazioneDStabilita = 0;    

% scegliere ampiezza impulso (dMAx = 95)
d              = 95;

%--------------------------------------------------------------------------
% Definizione Parametri Modello
M  = 3;                                       % Massa Cart 
m  = 10;                                      % Massa Pendolo al CG 
l  = 0.22;                                    % Altezza CG dal Cart 
g  = 9.81;                                    % Gravità
Wp = 0.20;                                    % Larghezza Piatto Pendolo
Hp = 0.03;                                    % Altezza Piatto Pendolo
Ip = (1/12)*(Wp^2 + Hp^2)*m;                  % Inerzia Rispetto a Yaw
mu = 0.90;                                    % Attrito Pneumatico - Strada

% Definizione Stato Iniziale Sistema
x0 = [0 0 0 0]';

% Parametri D-stabilità (Fissati Secondo Le Specifiche)
alfa  = 1.5;
theta = (53.84)*(pi/180);
raggio= 50;

% Parametri Simulazione
ts = 5;                                       % Tempo Di Simulazione
Ts = 0.01;                                    % Tempo Di Campionamento

%--------------------------------------------------------------------------
% Parametri Controllo Robusto
lmin = 0.005;                                 % valore minimo distanza CG cart
lmax = 0.35;                                  % valore massimo distanza CG cart
Ipmin= 0.005;                                 % valore minimo inerzia
Ipmax= 0.35;                                  % valore massimo inerzia

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

%--------------------------------------------------------------------------
% Determinazione Guadagni Sistema
switch controllerType
    case 1
        K = controlloreStabilizzante(A,B1,4,1);
    case 2
        K = controlloreDstabilita(A,B1,4,1,alfa,theta,raggio);
    case 3
        K = controlloreHinfinito(A,B1,B2,[1 0 0 0],0,0,4,1,abilitazioneDStabilita,alfa,theta,raggio);
    case 4
        K = controlloreH2(A,B1,B2,[1 0 0 0],0,4,1,abilitazioneDStabilita,alfa,theta,raggio);
    case 5
        K = controlloreL1(A,B1,B2,[1 0 0 0],0,0,4,1,abilitazioneDStabilita,alfa,theta,raggio);
    case 6
        KHinf = controlloreHinfinito(A,B1,B2,[1 0 0 0],0,0,4,1,abilitazioneDStabilita,alfa,theta,raggio);
        KH2   = controlloreH2(A,B1,B2,[1 0 0 0],0,4,1,abilitazioneDStabilita,alfa,theta,raggio);
        KL1   = controlloreL1(A,B1,B2,[1 0 0 0],0,0,4,1,abilitazioneDStabilita,alfa,theta,raggio);
    case 7
        [Ak,Bk,Ck,Dk] = determinaControlloreHinf(A,B1,B2,C,D,[1 0 0 0],0,0,4,1,2);
    case 8
        [Ak,Bk,Ck,Dk] = determinazioneControlloreH2(A,B1,B2,C,D,[1 0 0 0],0,4,1,2);
    case 9
        [Ak1,Bk1,Ck1,Dk1] = determinaControlloreHinf(A,B1,B2,C,D,[1 0 0 0],0,0,4,1,2);
        [Ak2,Bk2,Ck2,Dk2] = determinazioneControlloreH2(A,B1,B2,C,D,[1 0 0 0],0,4,1,2);
    case 10
        K = controlloreRobustoStabilizzante(A_f,B1_f,lmin,lmax,Ipmin,Ipmax,4,1);
    case 11
        K = controlloreRobustoHinf(A_f,B1_f,B2_f,lmin,lmax,Ipmin,Ipmax,[1 0 0 0],0,0,4,1);
    case 12
        K = controlloreRobustoH2(A_f,B1_f,B2_f,lmin,lmax,Ipmin,Ipmax,[1 0 0 0],0,4,1);
    case 13
        KinfR = controlloreRobustoHinf(A_f,B1_f,B2_f,lmin,lmax,Ipmin,Ipmax,[1 0 0 0],0,0,4,1);
        K2R   = controlloreRobustoH2(A_f,B1_f,B2_f,lmin,lmax,Ipmin,Ipmax,[1 0 0 0],0,4,1);
    case 14
        Kreg  = paretoOptimum(A,B1,B2,C,D,[1 0 0 0],[0 0 1 0],4,1);
    otherwise
        p     = [-2 -3 -4 -5];
        K     = -place(A,B1,p);
end

%--------------------------------------------------------------------------
% Rappresentazione Scelta Controllore
controllori = {"Stabilizzazione","D-Stabilità","Hinf Statico","H2 Statico","L1 Statico","Statici e Confronto","Hinf Dinamico","H2 Dinamico","Dinamici e Confronto","Stabilizzante Robusto","Hinf Robusto","H2 Robusto","Robusti e Confronto","Pareto"};

if controllerType >= 1 && controllerType <= 14
    string = strcat('-> Controllore Selezionato : ',controllori(controllerType));
    display(string);
else
    display("-> Controllore Selezionato : Pole Placement");
end

if controllerType == 3 || controllerType == 4 || controllerType == 5 
    if(abilitazioneDStabilita)
       display("-> D-Stabilità : On");
    else
       display("-> D-Stabilità : Off");
    end
end

%--------------------------------------------------------------------------
% Simulazione Sistema Reale
display("-> Simulazione Su Robot");
%--------------------------------------------------------------------------
% Progettazione stimatore stato (Luemberger)
p = [-10 -11 -12 -13];
H = place(A',C',p)';
As= A - H*C;
Bs=[B1 H];
Cs=[1 0 0 0;0 1 0 0;0 0 1 0;0 0 0 1];
Ds=[0 0 0;0 0 0;0 0 0;0 0 0];
%--------------------------------------------------------------------------
% Visualizzazione
ts          = 25;
t           = 0:Ts:ts;                              % intervallo simulazione
thetaSim    = zeros(length(t),1);                   % conservazione simulazione
thetaSimCom = zeros(length(t),3);            
% Simulazione Su Modello Non Lineare
switch controllerType
    case {1,2,3,4,5}
        sim("../SchemiSimulink/MyProject.slx");
        thetaSim(:,1) = ans.theta;
        figure(1);
        plot(t,thetaSim(:,1),'b','LineWidth',1);
        grid;
    case 6
        K = KHinf;
        sim("../SchemiSimulink/MyProject.slx");
        thetaSimCom(:,1) = ans.theta;
        K = KH2;
        sim("../SchemiSimulink/MyProject.slx");
        thetaSimCom(:,2) = ans.theta;
        K = KL1;
        sim("../SchemiSimulink/MyProject.slx");
        thetaSimCom(:,3) = ans.theta;
        figure(1);
        plot(t,thetaSimCom(:,1),'b',t,thetaSimCom(:,2),'r',t,thetaSimCom(:,3),'m','LineWidth',1);
        legend("Hinf","H2","L1");
        grid;
    case 7
        sim("../SchemiSimulink/MyProject1.slx");
        thetaSim(:,1) = ans.theta;
        figure(1);
        plot(t,thetaSim(:,1),'b','LineWidth',1);
        grid;
    case 8
        sim("../SchemiSimulink/MyProject1.slx");
        thetaSim(:,1) = ans.theta;
        figure(1);
        plot(t,thetaSim(:,1),'b','LineWidth',1);
        grid;
    case 9
        Ak = Ak1; Bk = Bk1; Ck = Ck1; Dk = Dk1;
        sim("../SchemiSimulink/MyProject1.slx");
        thetaSimCom(:,1) = ans.theta;
        Ak = Ak2; Bk = Bk2; Ck = Ck2; Dk = Dk2;
        sim("../SchemiSimulink/MyProject1.slx");
        thetaSimCom(:,2) = ans.theta;
        thetaSimCom(:,3) = ans.theta;
        figure(1);
        plot(t,thetaSimCom(:,1),'b',t,thetaSimCom(:,2),'r',t,thetaSimCom(:,3),'r','LineWidth',1);
        legend("Hinf","H2");
        grid;
    case 10
        sim("../SchemiSimulink/MyProject.slx");
        thetaSim(:,1) = ans.theta;
        figure(1);
        plot(t,thetaSim(:,1),'b','LineWidth',1);
        grid;
    case 11
        sim("../SchemiSimulink/MyProject.slx");
        thetaSim(:,1) = ans.theta;
        figure(1);
        plot(t,thetaSim(:,1),'b','LineWidth',1);
        grid;
    case 12
        sim("../SchemiSimulink/MyProject.slx");
        thetaSim(:,1) = ans.theta;
        figure(1);
        plot(t,thetaSim(:,1),'b','LineWidth',1);
        grid;
    case 13
        K = KinfR;
        sim("../SchemiSimulink/MyProject.slx");
        thetaSimCom(:,1) = ans.theta;
        K = K2R;
        sim("../SchemiSimulink/MyProject.slx");
        thetaSimCom(:,2) = ans.theta;
        thetaSimCom(:,3) = ans.theta;
        figure(1);
        plot(t,thetaSimCom(:,1),'b',t,thetaSimCom(:,2),'r',t,thetaSimCom(:,3),'r','LineWidth',1);
        legend("Hinf","H2");
        grid;
    case 14
        K = Kreg;
        thetaSim = zeros(length(t),2);
        sim("../SchemiSimulink/MyProject.slx");
        thetaSim(:,1) = ans.theta;
        thetaSim(:,2) = ans.xc;
        figure(1);
        plot(t,thetaSim(:,1),'b','LineWidth',1);
        legend("theta");
        grid;
        figure(2);
        plot(t,thetaSim(:,2),'b','LineWidth',1);
        legend("xc");
        grid;
    otherwise 
        sim("../SchemiSimulink/MyProject.slx");
        thetaSim(:,1) = ans.theta;
        figure(1);
        plot(t,thetaSim(:,1),'b','LineWidth',1);
        grid;
end
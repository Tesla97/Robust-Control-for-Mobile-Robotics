%% Matlab Code For An Adaptive Cruise Conroller - Nicola Corea 235279
clear all
close all
clc;
%--------------------------------------------------------------------------
% Parameter of The Model
Ts  = 0.1;
tau = 0.5;
%--------------------------------------------------------------------------
% State Space Model
A   = [0 1 0;0 0 1;0 0 (-1/tau)];
B   = [0 0 (1/tau)];
C   = [0 1 0];
D   = 0;
%--------------------------------------------------------------------------
% Control Parameters
h      = 1;
lambda = 50;
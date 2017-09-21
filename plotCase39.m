clear all; 
close all;
clc;
CurrentDirectory=pwd;
ResultsPath='Results/case39wmac_con/';
SsControlOptions={'LQR-OPF','ALQR-OPF','OPF'};
LfControlOptions={'LQR'};


cd(ResultsPath);


cd('LQR-OPF/LQR/'); 
LqrOpf=load('case39wmac_con_LQR-OPF_alphapoint6_LQR.mat'); 
cd('..'); 
cd('..'); 

cd('ALQR-OPF/LQR/'); 
ALqrOpf=load('case39wmac_con_ALQR-OPF_alphapoint6_LQR.mat'); 
cd('..'); 
cd('..'); 


cd('OPF/LQR/'); 
Opf=load('case39wmac_con_OPF_alphapoint6_LQR.mat'); 
cd('..'); 
cd('..');

t=LqrOpf.t;
TFinal=LqrOpf.TFinal;
G=LqrOpf.G;
OMEGAS=LqrOpf.OMEGAS;

cd(CurrentDirectory);


%% 






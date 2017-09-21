clear all; 
close all;
clc;
CurrentDirectory=pwd;
ResultsPath='Results/case39wmac_con/';
SsControlOptions={'LQR-OPF','ALQR-OPF','OPF'};
LfControlOptions={'AGC'};


cd(ResultsPath);


cd('LQR-OPF/AGC/'); 
LqrOpf=load('case39wmac_con_LQR-OPF_alphapoint6_AGC.mat'); 
cd('..'); 
cd('..'); 

cd('ALQR-OPF/AGC/'); 
ALqrOpf=load('case39wmac_con_ALQR-OPF_alphapoint6_AGC.mat'); 
cd('..'); 
cd('..'); 


cd('OPF/AGC/'); 
Opf=load('case39wmac_con_OPF_alphapoint6_AGC.mat'); 
cd('..'); 
cd('..');

t=LqrOpf.t;
TFinal=LqrOpf.TFinal;
G=LqrOpf.G;
OMEGAS=LqrOpf.OMEGAS;
AreaSet=LqrOpf.AreaSet;

cd(CurrentDirectory);


%% 






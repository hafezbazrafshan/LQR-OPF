clear all; 
close all;
clc;

CurrentDirectory=pwd;
% importing case200 LQR-OPF
cd('Results/case_illinois200/LQR-OPF/LQR');
LQROPF=load('case_illinois200_LQR-OPF_LQRalphapoint8.mat'); 
cd(CurrentDirectory); 

% importing case200 ALQR-OPF
cd('Results/case_illinois200/ALQR-OPF/LQR');
ALQROPF=load('case_illinois200_ALQR-OPF_LQRalphapoint8.mat'); 
cd(CurrentDirectory); 


% importing case200 OPF
cd('Results/case_illinois200/OPF/LQR');
OPF=load('case_illinois200_OPF_LQRalphapoint8.mat'); 
cd(CurrentDirectory); 
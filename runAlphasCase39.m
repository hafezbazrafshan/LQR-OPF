clear all;
clc;
CaseFiles={'case39wmac_con'};
% Perturbation
PRatio=0.1; 
QRatio=0.0484;



% Coupling parameter
AlphaVec=[0; 0.2; 0.4; 0.6;0.8];

SsControlOptions={'LQR-OPF','OPF'};
SteadyStateOutput=cell(length(AlphaVec),length(SsControlOptions));

% LQR time importance
Tlqr=1000;


%% Steady-state runs
if exist('Results')~=7
mkdir('Results');
end
SaveName=['SteadyStateAlphaComparison','Percent.txt'];
FileID=fopen(['Results/',SaveName],'w'); 
fprintf(FileID,'%-15s & %-15s & %-15s & %-15s & %-15s & %-15s & %-15s & %-15s\n',...
    'Network', 'SsMethod','Alpha','ObjValue', 'SsCost', 'StCostEst.','TotalCostEst.', 'CompTime');
for mm=1:length(AlphaVec)
 Alpha=AlphaVec(mm);

    CaseFile=CaseFiles{1};
for ii=1:length(SsControlOptions)
    SteadyStateCase=steadyStateDriver(CaseFile,SsControlOptions{ii},Alpha, Tlqr, PRatio, QRatio);  
    SteadyStateOutput{mm,ii}=SteadyStateCase;
     fprintf(FileID, '%-15s & %-15s & %-15.2f & %-15.2f & %-15.2f & %-15.2f & %-15.2f & %-15.2f  \n', ...
    CaseFile, SsControlOptions{ii},Alpha, SteadyStateCase.SsObjEst, SteadyStateCase.SsCost,...
   SteadyStateCase.TrCostEstimate, SteadyStateCase.TotalCostEstimate, SteadyStateCase.CompTime);
end
end
fclose(FileID);


%% Dynamical simulations
LfControlOptions={'LQR'};
DynamicOutput=cell(length(AlphaVec), length(SsControlOptions), length(LfControlOptions));
SaveName=['DynamicalAlphaComparison'    ,'Percent.txt'];
FileID=fopen(['Results/',SaveName],'w'); 
        fprintf(FileID,'%-15s & %-15s & %-15s & & %-15s & %-15s & %-15s & %-15s & %-15s \n',...
    'Network', 'SsMethod','Alpha', 'SsCost','StCost','TotalCost.', 'MaxFreqDev', 'MaxVoltDev');
for mm=1:length(AlphaVec)
    CaseFile=CaseFiles{1};
              LfControl=LfControlOptions{1};
    for ii=1:length(SsControlOptions)
    DynamicCase=dynamicDriver([SteadyStateOutput{mm,ii}.SteadyStatePath, '/',SteadyStateOutput{mm,ii}.SteadyStateFileName], LfControl,'YesPlots');
    DynamicOutput{mm,ii,1}=DynamicCase;
     fprintf(FileID, '%-15s & %-15s &%-15.2f & %-15.2f & %-15.2f & %-15.2f & %-15.4f & %-15.4f  \n', ...
    CaseFile, SsControlOptions{ii},DynamicCase.Alpha, DynamicCase.SsCost, DynamicCase.TrCost, DynamicCase.TotalCost,...
    DynamicCase.MaxFreqDev, DynamicCase.MaxVoltDev );
end
end
fclose(FileID);







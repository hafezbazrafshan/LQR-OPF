clear all;
clc;
casefile='case9wmac_con';

% Coupling parameter
Alpha=0.8;

% LQR time importance
Tlqr=1000;

% LfControl Method:
LfControl='AGC';

% StControlOptions={'LQR-OPF', 'ALQR-OPF','OPF'};
% StControlOptions={'OPF','LQR-OPF'};
StControlOptions={'OPF'};

Output=cell(size(StControlOptions));


if exist('Results')~=7
mkdir('Results');
end
% SaveName=['Case9ReportAGC',num2str(Alpha*100),'Percent.txt'];
% FileID=fopen(['Results/',SaveName],'w'); 
% fprintf(FileID,'%-15s & %-15s & %-15s & %-15s & %-15s & %-15s  & %-15s & %-15s & %-15s & %-15s \n',...
%     'Network', 'Method', 'SsObjEst.', 'SsCost.', 'StCostEst.', 'StCost.', 'TotCost.', 'CompTime', 'MaxFreqDev.', 'MaxVoltDev.');



for ii=1:length(StControlOptions)
    Output{ii}=workflow(casefile,StControlOptions{ii},LfControl,Alpha, 'WithPlots');  
%      fprintf(FileID, '%-15s & %-15s & %-15.2f & %-15.2f & %-15.2f & %-15.2f  & %-15.2f & %-15.2f & %-15.4f %-15.4f \n', ...
%     casefile, StControlOptions{ii}, Output{ii}.SsObjEst, Output{ii}.SsCost, Output{ii}.TrCostEstimate, Output{ii}.TrCost, Output{ii}.TotalCost, Output{ii}.CompTime, ...
%     max(max(abs(Output{ii}.omegaVec-Output{ii}.OMEGAS)))./(2*pi), max(max(abs(Output{ii}.vVec-repmat(Output{ii}.vS, 1,Output{ii}.NSamples)))));
end

fclose(FileID);
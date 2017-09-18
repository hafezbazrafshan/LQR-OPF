clear all;
clc;
casefile='case2869pegase';

% Coupling parameter
Alpha=0.8;

% LQR time importance
Tlqr=1000;

% LfControl Method:
% LfControl='AGC';
LfControl='LQR';

% StControlOptions={'LQR-OPF', 'ALQR-OPF','OPF'};
StControlOptions={'OPF','ALQR-OPF'};
% StControlOptions={'OPF'};
% StControlOptions={'LQR-OPF'};

Output=cell(size(StControlOptions));


if exist('Results')~=7
mkdir('Results');
end
SaveName=['Case2869Report',num2str(Alpha*100),'Percent.txt'];
FileID=fopen(['Results/',SaveName],'w'); 
fprintf(FileID,'%-15s & %-15s & %-15s & %-15s & %-15s & %-15s  & %-15s \n',...
    'Network', 'Method', 'SsObjEst.', 'SsCost.', 'StCost.', 'TotCost.', 'CompTime');



for ii=1:length(StControlOptions)
    Output{ii}=workflow(casefile,StControlOptions{ii},LfControl,Alpha, 'WithPlots');  
     fprintf(FileID, '%-15s & %-15s & %-15.2f & %-15.2f & %-15.2f  & %-15.2f & %-15.2f \n', ...
    casefile, StControlOptions{ii}, Output{ii}.SsObjEst, Output{ii}.SsCost, Output{ii}.TrCost, Output{ii}.TotalCost, Output{ii}.CompTime);
end

fclose(FileID);
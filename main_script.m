clear all;
clc;

%% Problem Parameters:

% Network names
% casefiles={'case9wmac_con'; 'case14wmac_con';'case39wmac_con','case57','case_illinois200'};
% casefiles={'case39wmac_con'};
% casefiles={'case_illinois200'};


% casefiles={'case57'};
% casefiles={'case9wmac_con'};

casefiles={'case9wmac_con'};
LfControl='AGC';
Alpha=0.8;


for case_index=1:length(casefiles)
    casefile=casefiles{case_index};
    lqropf=workflow(casefile,'ALQR-OPF',LfControl,Alpha, 'WithPlots');
        if lqropf.N<250
            fprintf(fileID, '%-20s & %-10s & %-10.2f & %-20.2f & %-20.2f  & %-20.4f  %-20.4f\n', ...
    casefile, 'lqr-opf', lqropf.SsCost, lqropf.TrCost, lqropf.SsCost+lqropf.TrCost, ...
    max(max(abs(lqropf.omegaVec-lqropf.OMEGAS)))./(2*pi), max(max(abs(lqropf.vVec-repmat(lqropf.vS, 1,lqropf.NSamples)))));
        else 
            fprintf(fileID, '%-20s & %-10s & %-10.2f & %-20.2f & %-20.2f  & %-20.4f  %-20.4f\n', ...
               casefile, 'lqr-opf', lqropf.SsCost, lqropf.TrCostEstimate, lqropf.SsCost+lqropf.TrCostEstimate, ...
   '----', '----'); 
        end
     opf=workflow(casefile,'OPF',LfControl,Alpha);
    if lqropf.N<250
            fprintf(fileID, '%-20s & %-10s & %-10.2f & %-20.2f & %-20.2f  & %-20.4f  %-20.4f\n', ...
    casefile, 'opf', opf.SsCost, opf.TrCost, opf.SsCost+opf.TrCost, ...
    max(max(abs(opf.omegaVec-opf.OMEGAS)))./(2*pi), max(max(abs(opf.vVec-repmat(opf.vS,1,opf.NSamples)))));
        else 
            fprintf(fileID, '%-20s & %-10s & %-10.2f & %-20.2f & %-20.2f  & %-20.4f  %-20.4f\n', ...
               casefile, 'opf', opf.SsCost, opf.TrCostEstimate, opf.SsCost+opf.TrCostEstimate, ...
   '----', '----'); 
        end
    
end

fclose(fileID);



CurrentFolder=pwd;
LQROPFAddress='C:\Users\aju084\Documents\MATLAB\LQR-OPF\Results\case39wmac_con\LQR-OPF\LQR';
cd(LQROPFAddress);
Alpha0=load('case39wmac_con_LQR-OPF_alphapoint0_LQR.mat'); 
Alpha2=load('case39wmac_con_LQR-OPF_alphapoint2_LQR.mat'); 
Alpha4=load('case39wmac_con_LQR-OPF_alphapoint4_LQR.mat'); 
Alpha6=load('case39wmac_con_LQR-OPF_alphapoint6_LQR.mat'); 
Alpha8=load('case39wmac_con_LQR-OPF_alphapoint8_LQR.mat');
cd(CurrentFolder);

OPFAddress='C:\Users\aju084\Documents\MATLAB\LQR-OPF\Results\case39wmac_con\OPF\LQR';
cd(OPFAddress);
Alpha0OPF=load('case39wmac_con_OPF_alphapoint0_LQR.mat'); 
Alpha2OPF=load('case39wmac_con_OPF_alphapoint2_LQR.mat'); 
Alpha4OPF=load('case39wmac_con_OPF_alphapoint4_LQR.mat'); 
Alpha6OPF=load('case39wmac_con_OPF_alphapoint6_LQR.mat'); 
Alpha8OPF=load('case39wmac_con_OPF_alphapoint8_LQR.mat');
cd(CurrentFolder);

AlphaVec=[0;0.2;0.4;0.6;0.8];

x0=0;
y0=1;
width=8;
height=5;
figure1=figure('Units','inches',...
'Position',[x0 y0 width height],...
'PaperPositionMode','auto');
set(figure1, 'Name', 'BarComparison');
LQROPFCoordinates=AlphaVec-0.03;
LQROPFMat = [Alpha0.SsCost Alpha0.TrCost; ...
    Alpha2.SsCost Alpha2.TrCost;
    Alpha4.SsCost Alpha4.TrCost; ...
    Alpha6.SsCost Alpha6.TrCost; ...
    Alpha8.SsCost Alpha8.TrCost; ];
LQROPFBarPlot=bar(LQROPFCoordinates,  LQROPFMat,'stacked');
set(LQROPFBarPlot,  'BarWidth',0.3);
LQROPFBarPlot(1).FaceColor=[0 0 0.5]; % blue
LQROPFBarPlot(2).FaceColor=	[100,206,210]/255; % cyan
 
set(gca,'nextplot','add') ;
OPFCoordinates=AlphaVec+0.03;
OpfMat = [Alpha0OPF.SsCost Alpha0OPF.TrCost; ...
    Alpha2OPF.SsCost Alpha2OPF.TrCost;
    Alpha4OPF.SsCost Alpha4OPF.TrCost; ...
    Alpha6OPF.SsCost Alpha6OPF.TrCost; ...
    Alpha8OPF.SsCost Alpha8OPF.TrCost;];

OPFBarPlot=bar(OPFCoordinates,   OpfMat,'stacked');
set(OPFBarPlot,  'BarWidth',0.3);
OPFBarPlot(1).FaceColor=[0 0 255]/255; % blue
OPFBarPlot(2).FaceColor=[135,206,250]/255; %light sky blue

fig_handle=gca;

fig_handle.XTick=AlphaVec;
fig_handle.YLim=[0  max( max( [sum(LQROPFMat,2); sum(OpfMat,2)]))+40000];
fig_handle.YTick=[0:10000:80000];
fig_handle.YTickLabel=[num2str(fig_handle.YTick.'/1000),repmat('k',length(fig_handle.YTick.'),1)];


set(gca,'nextplot','add') ;
hold on
linePlot=plot(AlphaVec,sum(OpfMat,2)-sum(LQROPFMat,2));
linePlot.LineStyle='-';
linePlot.LineWidth=3;
linePlot.Color=[1 0.7 0];
linePlot.Marker='s';
linePlot.MarkerSize=10;
linePlot.MarkerFaceColor='w';


set(fig_handle,'box','on');
set(fig_handle,'fontSize',20); 
set(fig_handle,'defaulttextinterpreter','latex');
xlabel('Coupling coefficient $\alpha$', 'FontWeight','bold');
 ylabel('Cost (\$)'); 
fig_handle.TickLabelInterpreter='latex';

yticks=fig_handle.YTick.';
legendText={'LQR-OPF steady-state costs'; 'LQR-OPF control costs';'OPF steady-state costs'; 'OPF control costs'; 'LQR-OPF savings'};
 lgd=legend(fig_handle,legendText);
 lgd.FontSize=12;
 lgd.Interpreter='Latex';
 lgd.Location='North';
 lgd.Orientation='Vertical';

 
  cd('Figures'); 
  print -dpdf test39LQROPFbarplot.pdf
print -depsc2 test39LQROPFbarplot
cd('..');
%% LQR-OPF Freq.
ACEMax=max(max([LqrOpf.ACEVec,ALqrOpf.ACEVec, Opf.ACEVec]))./(2*pi);
ACEMin=min(min([LqrOpf.ACEVec,ALqrOpf.ACEVec, Opf.ACEVec]))./(2*pi);
ACEOffSet=0.05;

x0=0;
y0=1;
width=8;
height=5;
figure1=figure('Units','inches',...
'Position',[x0 y0 width height],...
'PaperPositionMode','auto');
set(figure1, 'Name', 'LQR-OPF:genfreq');
plot(t,LqrOpf.ACEVec.'/(2*pi),'lineWidth',2);
FigHandle=gca; 
set(FigHandle,'box','on');
set(FigHandle,'fontSize',20); 
set(FigHandle,'defaulttextinterpreter','latex');
FigHandle.TickLabelInterpreter='latex';
ytickformat(FigHandle,'%.2f'); 
axis(FigHandle,[0 TFinal ACEMin-ACEOffSet ACEMax+ACEOffSet+0.1]);
set(FigHandle,'YTick',floor(ACEMin*100./5)*0.05:0.1:ceil(ACEMax*100./5)*0.05);
xlabel('Time (sec)', 'FontWeight','bold');
 ylabel('$\mathrm{ACE}$ (pu)'); 
 title('ACE per area');

 LegendText=cellstr([repmat('Area ', length(AreaSet),1),num2str([1:length(AreaSet)].')]);
 Lgd=legend(FigHandle,LegendText);
 Lgd.FontSize=12;
 Lgd.Interpreter='Latex';
 Lgd.Location='East';
 grid on;
SubFigHandle=axes('position',[0.3 0.3 0.3 .3]);
ZoomIndex = (t<=5) & (t>=0);
plot(SubFigHandle,t(ZoomIndex),LqrOpf.ACEVec(:,ZoomIndex).'/(2*pi)) % plot on new axes
axis(SubFigHandle,'tight'); 
SubFigHandle.XTick= 0: 1: 5;
set(SubFigHandle,'XGrid','on'); 
set(SubFigHandle,'YGrid','on'); 
ytickformat(SubFigHandle,'%.2f'); 


set(SubFigHandle,'box','on');
set(SubFigHandle,'fontSize',14); 
set(SubFigHandle,'defaulttextinterpreter','latex');
SubFigHandle.TickLabelInterpreter='latex';
 SubFigHandle.XLabel.String='Time (sec)';
 SubFigHandle.YLabel.String='$\mathrm{ACE}$ (pu)'; 
 x=[0.22 0.15];
 y=[0.55 0.60];
 annotation('textarrow',x,y)

 
 if exist('Figures')~=7
    mkdir('Figures'); 
 end

 cd('Figures'); 
  print -dpdf Case39_LQROPF_alphapoint6_AGC_ACE.pdf
print -depsc2 Case39_LQROPF_alphapoint6_AGC_ACE
cd('..');



%% ALQR-OPF Freq.


x0=0;
y0=1;
width=8;
height=5;
figure1=figure('Units','inches',...
'Position',[x0 y0 width height],...
'PaperPositionMode','auto');
set(figure1, 'Name', 'LQR-OPF:genfreq');
plot(t,ALqrOpf.ACEVec.'/(2*pi),'lineWidth',2);
FigHandle=gca; 
set(FigHandle,'box','on');
set(FigHandle,'fontSize',20); 
set(FigHandle,'defaulttextinterpreter','latex');
xlabel('Time (sec)', 'FontWeight','bold');
 ylabel('$\mathrm{ACE}$ (pu)'); 
  title('ACE per area');
FigHandle.TickLabelInterpreter='latex';
axis(FigHandle,[0 TFinal ACEMin-ACEOffSet ACEMax+ACEOffSet+0.1]);
set(FigHandle,'YTick',floor(ACEMin*100./5)*0.05:0.1:ceil(ACEMax*100./5)*0.05);
ytickformat(FigHandle,'%.2f'); 
 LegendText=cellstr([repmat('Area ', length(AreaSet),1),num2str([1:length(AreaSet)].')]);
 Lgd=legend(FigHandle,LegendText);
 Lgd.FontSize=12;
 Lgd.Interpreter='Latex';
 Lgd.Location='East';
 grid on;
SubFigHandle=axes('position',[0.3 0.3 0.3 .3]);
ZoomIndex = (t<=5) & (t>=0);
plot(SubFigHandle,t(ZoomIndex),ALqrOpf.ACEVec(:,ZoomIndex).'/(2*pi)) % plot on new axes
axis(SubFigHandle,'tight'); 
SubFigHandle.XTick= 0: 1: 5;
set(SubFigHandle,'XGrid','on'); 
set(SubFigHandle,'YGrid','on'); 
ytickformat(SubFigHandle,'%.2f'); 

set(SubFigHandle,'box','on');
set(SubFigHandle,'fontSize',14); 
set(SubFigHandle,'defaulttextinterpreter','latex');
SubFigHandle.TickLabelInterpreter='latex';
 SubFigHandle.XLabel.String='Time (sec)';
 SubFigHandle.YLabel.String='$\mathrm{ACE}$ (pu)'; 
 x=[0.22 0.15];
 y=[0.55 0.60];
 annotation('textarrow',x,y)

 
 if exist('Figures')~=7
    mkdir('Figures'); 
 end

 cd('Figures'); 
  print -dpdf Case39_ALQROPF_alphapoint6_AGC_ACE.pdf
print -depsc2 Case39_ALQROPF_alphapoint6_AGC_ACE
cd('..');


%%

x0=0;
y0=1;
width=8;
height=5;
figure1=figure('Units','inches',...
'Position',[x0 y0 width height],...
'PaperPositionMode','auto');
set(figure1, 'Name', 'LQR-OPF:genfreq');
plot(t,Opf.ACEVec.'/(2*pi),'lineWidth',2);
FigHandle=gca; 
set(FigHandle,'box','on');
set(FigHandle,'fontSize',20); 
set(FigHandle,'defaulttextinterpreter','latex');
xlabel('Time (sec)', 'FontWeight','bold');
 ylabel('$\mathrm{ACE}$ (pu)');
  title('ACE per area');
FigHandle.TickLabelInterpreter='latex';
axis(FigHandle,[0 TFinal ACEMin-ACEOffSet ACEMax+ACEOffSet+0.1]);
set(FigHandle,'YTick',floor(ACEMin*100./5)*0.05:0.1:ceil(ACEMax*100./5)*0.05);
ytickformat(FigHandle,'%.2f'); 
 LegendText=cellstr([repmat('Area ', length(AreaSet),1),num2str([1:length(AreaSet)].')]);
 Lgd=legend(FigHandle,LegendText);
 Lgd.FontSize=12;
 Lgd.Interpreter='Latex';
 Lgd.Location='East';
 grid on;
SubFigHandle=axes('position',[0.3 0.3 0.3 .3]);
ZoomIndex = (t<=5) & (t>=0);
plot(SubFigHandle,t(ZoomIndex),Opf.ACEVec(:,ZoomIndex).'/(2*pi)) % plot on new axes
axis(SubFigHandle,'tight'); 
SubFigHandle.XTick= 0: 1: 5;
set(SubFigHandle,'XGrid','on'); 
set(SubFigHandle,'YGrid','on'); 
ytickformat(SubFigHandle,'%.2f'); 

set(SubFigHandle,'box','on');
set(SubFigHandle,'fontSize',14); 
set(SubFigHandle,'defaulttextinterpreter','latex');
SubFigHandle.TickLabelInterpreter='latex';
 SubFigHandle.XLabel.String='Time (sec)';
 SubFigHandle.YLabel.String='$\mathrm{ACE}$ (pu)'; 
 x=[0.22 0.15];
 y=[0.55 0.60];
 annotation('textarrow',x,y)

 
 if exist('Figures')~=7
    mkdir('Figures'); 
 end

 cd('Figures'); 
  print -dpdf Case39_OPF_alphapoint6_AGC_ACE.pdf
print -depsc2 Case39_OPF_alphapoint6_AGC_ACE
cd('..');
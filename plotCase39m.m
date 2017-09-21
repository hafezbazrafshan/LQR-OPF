%%  mech plots
MyMax=max(max([LqrOpf.mVec-repmat(LqrOpf.mS,1,length(t)),ALqrOpf.mVec-repmat(ALqrOpf.mS,1,length(t)),...
    Opf.mVec-repmat(Opf.mS,1,length(t))]));
MyMin=min(min([LqrOpf.mVec-repmat(LqrOpf.mS,1,length(t)),ALqrOpf.mVec-repmat(ALqrOpf.mS,1,length(t)),...
    Opf.mVec-repmat(Opf.mS,1,length(t))]));

MyOffSet=1;

x0=0;
y0=1;
width=8;
height=5;
Figure4=figure('Units','inches',...
'Position',[x0 y0 width height],...
'PaperPositionMode','auto');
set(Figure4, 'Name', 'LQR-OPF:mech');
plot(t,(LqrOpf.mVec-repmat(LqrOpf.mS,1,length(t))).','lineWidth',2);
FigHandle=gca; 
set(FigHandle,'box','on');
set(FigHandle,'fontSize',20); 
set(FigHandle,'defaulttextinterpreter','latex');
FigHandle.TickLabelInterpreter='latex';
ytickformat(FigHandle,'%.2f'); 
grid on;
xlabel('Time (sec)', 'FontWeight','bold');
 ylabel('$\mathbf{m}-\mathbf{m}^\mathrm{eq}$ (pu)'); 
  title('Generator mech. input');

axis(FigHandle,[0 TFinal MyMin-MyOffSet MyMax+MyOffSet]);
set(FigHandle,'YTick',floor(MyMin):1:ceil(MyMax));
 LegendText=cellstr([repmat('Gen. ', G,1),num2str([1:G].')]);
 Lgd=legend(FigHandle,LegendText);
 Lgd.FontSize=12;
 Lgd.Interpreter='Latex';
 Lgd.Location='East';
 grid on;


  SubFigHandle=axes('position',[0.3 0.3 0.3 .3]);
ZoomIndex = (t<=5) & (t>=0);
plot(SubFigHandle,t(ZoomIndex),(LqrOpf.mVec(:,ZoomIndex)-repmat(LqrOpf.mS,1,length(t(ZoomIndex)))).') % plot on new axes
axis(SubFigHandle,'tight'); 
SubFigHandle.XTick= 0: 1: 5;
ytickformat(SubFigHandle,'%.2f');
set(SubFigHandle,'box','on');
set(SubFigHandle,'fontSize',14); 
set(SubFigHandle,'defaulttextinterpreter','latex');
SubFigHandle.TickLabelInterpreter='latex';
 SubFigHandle.XLabel.String='Time (sec)';
 SubFigHandle.YLabel.String='$\mathbf{m}-\mathbf{m}^\mathrm{eq}$ (pu)'; 
 set(SubFigHandle,'XGrid','on'); 
set(SubFigHandle,'YGrid','on'); 
 x=[0.22 0.15];
 y=[0.55 0.60];
 annotation('textarrow',x,y)

 

 
 if exist('Figures')~=7
    mkdir('Figures'); 
 end

 cd('Figures'); 
  print -dpdf Case39_LQROPF_alphapoint6_LQR_m.pdf
print -depsc2 Case39_LQROPF_alphapoint6_LQR_m
cd('..');



%% ALQR-OPF

x0=0;
y0=1;
width=8;
height=5;
Figure4=figure('Units','inches',...
'Position',[x0 y0 width height],...
'PaperPositionMode','auto');
set(Figure4, 'Name', 'LQR-OPF:mech');
plot(t,(ALqrOpf.mVec-repmat(ALqrOpf.mS,1,length(t))).','lineWidth',2);
FigHandle=gca; 
set(FigHandle,'box','on');
set(FigHandle,'fontSize',20); 
set(FigHandle,'defaulttextinterpreter','latex');
FigHandle.TickLabelInterpreter='latex';
ytickformat(FigHandle,'%.2f'); 
grid on;
xlabel('Time (sec)', 'FontWeight','bold');
 ylabel('$\mathbf{m}-\mathbf{m}^\mathrm{eq}$ (pu)'); 
  title('Generator mech. input');

FigHandle.TickLabelInterpreter='latex';
axis(FigHandle,[0 TFinal MyMin-MyOffSet MyMax+MyOffSet]);
set(FigHandle,'YTick',floor(MyMin):1:ceil(MyMax));
 LegendText=cellstr([repmat('Gen. ', G,1),num2str([1:G].')]);
 Lgd=legend(FigHandle,LegendText);
 Lgd.FontSize=12;
 Lgd.Interpreter='Latex';
 Lgd.Location='East';
 grid on;


  SubFigHandle=axes('position',[0.3 0.3 0.3 .3]);
ZoomIndex = (t<=5) & (t>=0);
plot(SubFigHandle,t(ZoomIndex),(ALqrOpf.mVec(:,ZoomIndex)-repmat(ALqrOpf.mS,1,length(t(ZoomIndex)))).') % plot on new axes
axis(SubFigHandle,'tight'); 
SubFigHandle.XTick= 0: 1: 5;
ytickformat(SubFigHandle,'%.2f');
set(SubFigHandle,'box','on');
set(SubFigHandle,'fontSize',14); 
set(SubFigHandle,'defaulttextinterpreter','latex');
SubFigHandle.TickLabelInterpreter='latex';
 SubFigHandle.XLabel.String='Time (sec)';
 SubFigHandle.YLabel.String='$\mathbf{m}-\mathbf{m}^\mathrm{eq}$ (pu)'; 
 set(SubFigHandle,'XGrid','on'); 
set(SubFigHandle,'YGrid','on'); 
 x=[0.22 0.15];
 y=[0.55 0.60];
 annotation('textarrow',x,y)

 

 
 if exist('Figures')~=7
    mkdir('Figures'); 
 end

 cd('Figures'); 
  print -dpdf Case39_ALQROPF_alphapoint6_LQR_m.pdf
print -depsc2 Case39_ALQROPF_alphapoint6_LQR_m
cd('..');

%% mech plots


x0=0;
y0=1;
width=8;
height=5;
Figure4=figure('Units','inches',...
'Position',[x0 y0 width height],...
'PaperPositionMode','auto');
set(Figure4, 'Name', 'LQR-OPF:mech');
plot(t,(Opf.mVec-repmat(Opf.mS,1,length(t))).','lineWidth',2);
FigHandle=gca; 
set(FigHandle,'box','on');
set(FigHandle,'fontSize',20); 
set(FigHandle,'defaulttextinterpreter','latex');
FigHandle.TickLabelInterpreter='latex';
ytickformat(FigHandle,'%.2f'); 
grid on;
xlabel('Time (sec)', 'FontWeight','bold');
 ylabel('$\mathbf{m}-\mathbf{m}^\mathrm{eq}$ (pu)'); 
  title('Generator mech. input');

FigHandle.TickLabelInterpreter='latex';
axis(FigHandle,[0 TFinal MyMin-MyOffSet MyMax+MyOffSet]);
set(FigHandle,'YTick',floor(MyMin):1:ceil(MyMax));
 LegendText=cellstr([repmat('Gen. ', G,1),num2str([1:G].')]);
 Lgd=legend(FigHandle,LegendText);
 Lgd.FontSize=12;
 Lgd.Interpreter='Latex';
 Lgd.Location='East';
 grid on;


  SubFigHandle=axes('position',[0.3 0.3 0.3 .3]);
ZoomIndex = (t<=5) & (t>=0);
plot(SubFigHandle,t(ZoomIndex),(Opf.mVec(:,ZoomIndex)-repmat(Opf.mS,1,length(t(ZoomIndex)))).') % plot on new axes
axis(SubFigHandle,'tight'); 
SubFigHandle.XTick= 0: 1: 5;
ytickformat(SubFigHandle,'%.2f'); 

set(SubFigHandle,'box','on');
set(SubFigHandle,'fontSize',14); 
set(SubFigHandle,'defaulttextinterpreter','latex');
SubFigHandle.TickLabelInterpreter='latex';
 SubFigHandle.XLabel.String='Time (sec)';
 SubFigHandle.YLabel.String='$\mathbf{m}-\mathbf{m}^\mathrm{eq}$ (pu)'; 
 set(SubFigHandle,'XGrid','on'); 
set(SubFigHandle,'YGrid','on'); 
 x=[0.22 0.15];
 y=[0.55 0.60];
 annotation('textarrow',x,y)

 

 
 if exist('Figures')~=7
    mkdir('Figures'); 
 end

 cd('Figures'); 
  print -dpdf Case39_OPF_alphapoint6_LQR_m.pdf
print -depsc2 Case39_OPF_alphapoint6_LQR_m
cd('..');




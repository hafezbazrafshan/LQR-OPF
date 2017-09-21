%%  mech plots
MyMax=max(max([LqrOpf.fVec-repmat(LqrOpf.fS,1,length(t)),ALqrOpf.fVec-repmat(ALqrOpf.fS,1,length(t)),...
    Opf.fVec-repmat(Opf.fS,1,length(t))]));
MyMin=min(min([LqrOpf.fVec-repmat(LqrOpf.fS,1,length(t)),ALqrOpf.fVec-repmat(ALqrOpf.fS,1,length(t)),...
    Opf.fVec-repmat(Opf.fS,1,length(t))]));

MyOffSet=1;

x0=0;
y0=1;
width=8;
height=5;
Figure4=figure('Units','inches',...
'Position',[x0 y0 width height],...
'PaperPositionMode','auto');
set(Figure4, 'Name', 'LQR-OPF:f');
plot(t,(LqrOpf.fVec-repmat(LqrOpf.fS,1,length(t))).','lineWidth',2);
FigHandle=gca; 
set(FigHandle,'box','on');
set(FigHandle,'fontSize',20); 
set(FigHandle,'defaulttextinterpreter','latex');
grid on;
xlabel('Time (sec)', 'FontWeight','bold');
 ylabel('$\mathbf{f}-\mathbf{f}^\mathrm{eq}$ (pu)'); 
  title('Exciter field voltage');

FigHandle.TickLabelInterpreter='latex';
axis(FigHandle,[0 TFinal MyMin-MyOffSet MyMax+MyOffSet]);
set(FigHandle,'YTick',floor(MyMin):1:ceil(MyMax));
 LegendText=cellstr([repmat('Gen. ', G,1),num2str([1:G].')]);
 Lgd=legend(FigHandle,LegendText);
 Lgd.FontSize=12;
 Lgd.Interpreter='Latex';
 Lgd.Location='East';
 grid on;


  SubFigHandle=axes('position',[0.3 0.57 0.3 .3]);
ZoomIndex = (t<=5) & (t>=0);
plot(SubFigHandle,t(ZoomIndex),(LqrOpf.fVec(:,ZoomIndex)-repmat(LqrOpf.fS,1,length(t(ZoomIndex)))).') % plot on new axes
axis(SubFigHandle,'tight'); 
SubFigHandle.XTick= 0: 1: 5;

set(SubFigHandle,'box','on');
set(SubFigHandle,'fontSize',14); 
set(SubFigHandle,'defaulttextinterpreter','latex');
SubFigHandle.TickLabelInterpreter='latex';
 SubFigHandle.XLabel.String='Time (sec)';
 SubFigHandle.YLabel.String='$\mathbf{f}-\mathbf{f}^\mathrm{eq}$ (pu)'; 
 x=[0.22 0.15];
 y=[0.55 0.50];
 annotation('textarrow',x,y)

 

 
 if exist('Figures')~=7
    mkdir('Figures'); 
 end

 cd('Figures'); 
  print -dpdf Case39_LQROPF_alphapoint6_LQR_f.pdf
print -depsc2 Case39_LQROPF_alphapoint6_LQR_f
cd('..');



%% ALQR-OPF

x0=0;
y0=1;
width=8;
height=5;
Figure4=figure('Units','inches',...
'Position',[x0 y0 width height],...
'PaperPositionMode','auto');
set(Figure4, 'Name', 'LQR-OPF:f');
plot(t,(ALqrOpf.fVec-repmat(ALqrOpf.fS,1,length(t))).','lineWidth',2);
FigHandle=gca; 
set(FigHandle,'box','on');
set(FigHandle,'fontSize',20); 
set(FigHandle,'defaulttextinterpreter','latex');
grid on;
xlabel('Time (sec)', 'FontWeight','bold');
 ylabel('$\mathbf{f}-\mathbf{f}^\mathrm{eq}$ (pu)'); 
  title('Exciter field voltage');

FigHandle.TickLabelInterpreter='latex';
axis(FigHandle,[0 TFinal MyMin-MyOffSet MyMax+MyOffSet]);
set(FigHandle,'YTick',floor(MyMin):1:ceil(MyMax));
 LegendText=cellstr([repmat('Gen. ', G,1),num2str([1:G].')]);
 Lgd=legend(FigHandle,LegendText);
 Lgd.FontSize=12;
 Lgd.Interpreter='Latex';
 Lgd.Location='East';
 grid on;


  SubFigHandle=axes('position',[0.3 0.57 0.3 .3]);
ZoomIndex = (t<=5) & (t>=0);
plot(SubFigHandle,t(ZoomIndex),(ALqrOpf.fVec(:,ZoomIndex)-repmat(ALqrOpf.fS,1,length(t(ZoomIndex)))).') % plot on new axes
axis(SubFigHandle,'tight'); 
SubFigHandle.XTick= 0: 1: 5;

set(SubFigHandle,'box','on');
set(SubFigHandle,'fontSize',14); 
set(SubFigHandle,'defaulttextinterpreter','latex');
SubFigHandle.TickLabelInterpreter='latex';
 SubFigHandle.XLabel.String='Time (sec)';
 SubFigHandle.YLabel.String='$\mathbf{f}-\mathbf{f}^\mathrm{eq}$ (pu)'; 
 x=[0.22 0.15];
 y=[0.55 0.50];
 annotation('textarrow',x,y)

 

 
 if exist('Figures')~=7
    mkdir('Figures'); 
 end

 cd('Figures'); 
  print -dpdf Case39_ALQROPF_alphapoint6_LQR_f.pdf
print -depsc2 Case39_ALQROPF_alphapoint6_LQR_f
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
plot(t,(Opf.fVec-repmat(Opf.fS,1,length(t))).','lineWidth',2);
FigHandle=gca; 
set(FigHandle,'box','on');
set(FigHandle,'fontSize',20); 
set(FigHandle,'defaulttextinterpreter','latex');
grid on;
xlabel('Time (sec)', 'FontWeight','bold');
 ylabel('$\mathbf{f}-\mathbf{f}^\mathrm{eq}$ (pu)'); 
  title('Exciter field voltage');

FigHandle.TickLabelInterpreter='latex';
axis(FigHandle,[0 TFinal MyMin-MyOffSet MyMax+MyOffSet]);
set(FigHandle,'YTick',floor(MyMin):1:ceil(MyMax));
 LegendText=cellstr([repmat('Gen. ', G,1),num2str([1:G].')]);
 Lgd=legend(FigHandle,LegendText);
 Lgd.FontSize=12;
 Lgd.Interpreter='Latex';
 Lgd.Location='East';
 grid on;


  SubFigHandle=axes('position',[0.3 0.57 0.3 .3]);
ZoomIndex = (t<=5) & (t>=0);
plot(SubFigHandle,t(ZoomIndex),(Opf.fVec(:,ZoomIndex)-repmat(Opf.fS,1,length(t(ZoomIndex)))).') % plot on new axes
axis(SubFigHandle,'tight'); 
SubFigHandle.XTick= 0: 1: 5;

set(SubFigHandle,'box','on');
set(SubFigHandle,'fontSize',14); 
set(SubFigHandle,'defaulttextinterpreter','latex');
SubFigHandle.TickLabelInterpreter='latex';
 SubFigHandle.XLabel.String='Time (sec)';
 SubFigHandle.YLabel.String='$\mathbf{f}-\mathbf{f}^\mathrm{eq}$ (pu)'; 
 x=[0.22 0.15];
 y=[0.55 0.50];
 annotation('textarrow',x,y)

 

 
 if exist('Figures')~=7
    mkdir('Figures'); 
 end

 cd('Figures'); 
  print -dpdf Case39_OPF_alphapoint6_LQR_f.pdf
print -depsc2 Case39_OPF_alphapoint6_LQR_f
cd('..');




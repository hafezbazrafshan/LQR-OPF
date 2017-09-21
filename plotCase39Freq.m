%% LQR-OPF Freq.
FreqyMax=max(max([LqrOpf.omegaVec,ALqrOpf.omegaVec, Opf.omegaVec]))./(2*pi);
FreqyMin=min(min([LqrOpf.omegaVec,ALqrOpf.omegaVec, Opf.omegaVec]))./(2*pi);
FreqyOffSet=0.05;

x0=0;
y0=1;
width=8;
height=5;
figure1=figure('Units','inches',...
'Position',[x0 y0 width height],...
'PaperPositionMode','auto');
set(figure1, 'Name', 'LQR-OPF:genfreq');
plot(t,LqrOpf.omegaVec.'/(2*pi),'lineWidth',2);
FigHandle=gca; 
set(FigHandle,'box','on');
set(FigHandle,'fontSize',20); 
set(FigHandle,'defaulttextinterpreter','latex');
xlabel('Time (sec)', 'FontWeight','bold');
 ylabel('$\frac{1}{2\pi}\bf{\omega}$ (Hz)'); 
 title('Generator frequencies');
FigHandle.TickLabelInterpreter='latex';
axis(FigHandle,[0 TFinal FreqyMin-FreqyOffSet FreqyMax+FreqyOffSet+0.1]);
set(FigHandle,'YTick',floor(FreqyMin*100./5)*0.05:0.05:ceil(FreqyMax*100./5)*0.05);
 LegendText=cellstr([repmat('Gen. ', G,1),num2str([1:G].')]);
 Lgd=legend(FigHandle,LegendText);
 Lgd.FontSize=12;
 Lgd.Interpreter='Latex';
 Lgd.Location='East';
 grid on;
SubFigHandle=axes('position',[0.26 0.6 0.3 .3]);
ZoomIndex = (t<=5) & (t>=0);
plot(SubFigHandle,t(ZoomIndex),LqrOpf.omegaVec(:,ZoomIndex).'/(2*pi)) % plot on new axes
axis(SubFigHandle,'tight'); 
SubFigHandle.XTick= 0: 1: 5;
set(SubFigHandle,'XGrid','on'); 
set(SubFigHandle,'YGrid','on'); 

set(SubFigHandle,'box','on');
set(SubFigHandle,'fontSize',14); 
set(SubFigHandle,'defaulttextinterpreter','latex');
SubFigHandle.TickLabelInterpreter='latex';
 SubFigHandle.XLabel.String='Time (sec)';
 SubFigHandle.YLabel.String='$\frac{1}{2\pi}\bf{\omega}$ (Hz)'; 
 x=[0.21 0.16];
 y=[0.55 0.49];
 annotation('textarrow',x,y)

 
 if exist('Figures')~=7
    mkdir('Figures'); 
 end

 cd('Figures'); 
  print -dpdf Case39_LQROPF_alphapoint6_LQR_freq.pdf
print -depsc2 Case39_LQROPF_alphapoint6_LQR_freq
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
plot(t,ALqrOpf.omegaVec.'/(2*pi),'lineWidth',2);
FigHandle=gca; 
set(FigHandle,'box','on');
set(FigHandle,'fontSize',20); 
set(FigHandle,'defaulttextinterpreter','latex');
xlabel('Time (sec)', 'FontWeight','bold');
 ylabel('$\frac{1}{2\pi}\bf{\omega}$ (Hz)'); 
  title('Generator frequencies');
FigHandle.TickLabelInterpreter='latex';
axis(FigHandle,[0 TFinal FreqyMin-FreqyOffSet FreqyMax+FreqyOffSet+0.1]);
set(FigHandle,'YTick',floor(FreqyMin*100./5)*0.05:0.05:ceil(FreqyMax*100./5)*0.05);
 LegendText=cellstr([repmat('Gen. ', G,1),num2str([1:G].')]);
 Lgd=legend(FigHandle,LegendText);
 Lgd.FontSize=12;
 Lgd.Interpreter='Latex';
 Lgd.Location='East';
 grid on;
SubFigHandle=axes('position',[0.26 0.6 0.3 .3]);
ZoomIndex = (t<=5) & (t>=0);
plot(SubFigHandle,t(ZoomIndex),ALqrOpf.omegaVec(:,ZoomIndex).'/(2*pi)) % plot on new axes
axis(SubFigHandle,'tight'); 
SubFigHandle.XTick= 0: 1: 5;
set(SubFigHandle,'XGrid','on'); 
set(SubFigHandle,'YGrid','on'); 

set(SubFigHandle,'box','on');
set(SubFigHandle,'fontSize',14); 
set(SubFigHandle,'defaulttextinterpreter','latex');
SubFigHandle.TickLabelInterpreter='latex';
 SubFigHandle.XLabel.String='Time (sec)';
 SubFigHandle.YLabel.String='$\frac{1}{2\pi}\bf{\omega}$ (Hz)'; 
 x=[0.21 0.16];
 y=[0.55 0.49];
 annotation('textarrow',x,y)

 
 if exist('Figures')~=7
    mkdir('Figures'); 
 end

 cd('Figures'); 
  print -dpdf Case39_ALQROPF_alphapoint6_LQR_freq.pdf
print -depsc2 Case39_ALQROPF_alphapoint6_LQR_freq
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
plot(t,Opf.omegaVec.'/(2*pi),'lineWidth',2);
FigHandle=gca; 
set(FigHandle,'box','on');
set(FigHandle,'fontSize',20); 
set(FigHandle,'defaulttextinterpreter','latex');
xlabel('Time (sec)', 'FontWeight','bold');
 ylabel('$\frac{1}{2\pi}\bf{\omega}$ (Hz)');
  title('Generator frequencies');
FigHandle.TickLabelInterpreter='latex';
axis(FigHandle,[0 TFinal FreqyMin-FreqyOffSet FreqyMax+FreqyOffSet+0.1]);
set(FigHandle,'YTick',floor(FreqyMin*100./5)*0.05:0.05:ceil(FreqyMax*100./5)*0.05);
 LegendText=cellstr([repmat('Gen. ', G,1),num2str([1:G].')]);
 Lgd=legend(FigHandle,LegendText);
 Lgd.FontSize=12;
 Lgd.Interpreter='Latex';
 Lgd.Location='East';
 grid on;
SubFigHandle=axes('position',[0.26 0.6 0.3 .3]);
ZoomIndex = (t<=5) & (t>=0);
plot(SubFigHandle,t(ZoomIndex),Opf.omegaVec(:,ZoomIndex).'/(2*pi)) % plot on new axes
axis(SubFigHandle,'tight'); 
SubFigHandle.XTick= 0: 1: 5;
set(SubFigHandle,'XGrid','on'); 
set(SubFigHandle,'YGrid','on'); 

set(SubFigHandle,'box','on');
set(SubFigHandle,'fontSize',14); 
set(SubFigHandle,'defaulttextinterpreter','latex');
SubFigHandle.TickLabelInterpreter='latex';
 SubFigHandle.XLabel.String='Time (sec)';
 SubFigHandle.YLabel.String='$\frac{1}{2\pi}\bf{\omega}$ (Hz)'; 
 x=[0.21 0.16];
 y=[0.55 0.49];
 annotation('textarrow',x,y)

 
 if exist('Figures')~=7
    mkdir('Figures'); 
 end

 cd('Figures'); 
  print -dpdf Case39_OPF_alphapoint6_LQR_freq.pdf
print -depsc2 Case39_OPF_alphapoint6_LQR_freq
cd('..');
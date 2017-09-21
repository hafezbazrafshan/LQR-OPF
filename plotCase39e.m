%% e plots
eyMax=max(max([LqrOpf.eVec-repmat(LqrOpf.eS,1,length(t)),...
    ALqrOpf.eVec-repmat(ALqrOpf.eS,1,length(t)), Opf.eVec-repmat(Opf.eS,1,length(t))]));
eyMin=min(min([LqrOpf.eVec-repmat(LqrOpf.eS,1,length(t)),...
    ALqrOpf.eVec-repmat(ALqrOpf.eS,1,length(t)), Opf.eVec-repmat(Opf.eS,1,length(t))]));

eyOffSet=0.01;

x0=0;
y0=1;
width=8;
height=5;
Figure2=figure('Units','inches',...
'Position',[x0 y0 width height],...
'PaperPositionMode','auto');
set(Figure2, 'Name', 'GenEMF');
plot(t,LqrOpf.eVec-repmat(LqrOpf.eS,1,length(t)),'lineWidth',2);
FigHandle=gca; 
set(FigHandle,'box','on');
set(FigHandle,'fontSize',20); 
set(FigHandle,'defaulttextinterpreter','latex');
FigHandle.TickLabelInterpreter='latex';
ytickformat(FigHandle,'%.2f'); 
 xlabel('Time (sec)', 'FontWeight','bold');
 ylabel('$\boldmath{e}-\boldmath{e}^s$ (pu)'); 
axis([0 TFinal eyMin-eyOffSet eyMax+eyOffSet]);
xlabel('Time (sec)', 'FontWeight','bold');
 grid on;
title('Generator EMF dev.'); 
if exist('Figures')~=7
    mkdir('Figures'); 
 end

 cd('Figures'); 
  print -dpdf Case39_LQROPF_alphapoint6_LQR_e.pdf
print -depsc2 Case39_LQROPF_alphapoint6_LQR_e
cd('..');



%% ALQR-OPF

x0=0;
y0=1;
width=8;
height=5;
Figure2=figure('Units','inches',...
'Position',[x0 y0 width height],...
'PaperPositionMode','auto');
set(Figure2, 'Name', 'GenEMF');
plot(t,ALqrOpf.eVec-repmat(ALqrOpf.eS,1,length(t)),'lineWidth',2);
FigHandle=gca; 
set(FigHandle,'box','on');
set(FigHandle,'fontSize',20); 
set(FigHandle,'defaulttextinterpreter','latex');
FigHandle.TickLabelInterpreter='latex';
ytickformat(FigHandle,'%.2f'); 
 xlabel('Time (sec)', 'FontWeight','bold');
 ylabel('$\boldmath{e}-\boldmath{e}^s$ (pu)'); 
axis([0 TFinal eyMin-eyOffSet eyMax+eyOffSet]);
 grid on;
title('Generator EMF dev.'); 
if exist('Figures')~=7
    mkdir('Figures'); 
 end

 cd('Figures'); 
  print -dpdf Case39_ALQROPF_alphapoint6_LQR_e.pdf
print -depsc2 Case39_ALQROPF_alphapoint6_LQR_e
cd('..');


 
%% OPF

x0=0;
y0=1;
width=8;
height=5;
Figure2=figure('Units','inches',...
'Position',[x0 y0 width height],...
'PaperPositionMode','auto');
set(Figure2, 'Name', 'GenEMF');
plot(t,Opf.eVec-repmat(Opf.eS,1,length(t)),'lineWidth',2);
FigHandle=gca; 
set(FigHandle,'box','on');
set(FigHandle,'fontSize',20); 
set(FigHandle,'defaulttextinterpreter','latex');
FigHandle.TickLabelInterpreter='latex';
ytickformat(FigHandle,'%.2f'); 
 xlabel('Time (sec)', 'FontWeight','bold');
 ylabel('$\boldmath{e}-\boldmath{e}^s$ (pu)'); 
axis([0 TFinal eyMin-eyOffSet eyMax+eyOffSet]);
 grid on;
title('Generator EMF dev.'); 
if exist('Figures')~=7
    mkdir('Figures'); 
 end

 cd('Figures'); 
  print -dpdf Case39_OPF_alphapoint6_LQR_e.pdf
print -depsc2 Case39_OPF_alphapoint6_LQR_e
cd('..');





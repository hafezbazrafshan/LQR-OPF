%%  mech plots
MyMax=max(max([LqrOpf.vVec,ALqrOpf.vVec,...
    Opf.vVec]));
MyMin=min(min([LqrOpf.vVec,ALqrOpf.vVec,...
    Opf.vVec]));

MyOffSet=0.02;

x0=0;
y0=1;
width=8;
height=5;
Figure4=figure('Units','inches',...
'Position',[x0 y0 width height],...
'PaperPositionMode','auto');
set(Figure4, 'Name', 'LQR-OPF:v');
plot(t,(LqrOpf.vVec).','lineWidth',2);
FigHandle=gca; 
set(FigHandle,'box','on');
set(FigHandle,'fontSize',20); 
set(FigHandle,'defaulttextinterpreter','latex');
ytickformat(FigHandle,'%.2f'); 

grid on;
xlabel('Time (sec)', 'FontWeight','bold');
 ylabel('$\mathbf{v}-\mathbf{v}^\mathrm{eq}$ (pu)'); 
  title('Nodal voltages');

FigHandle.TickLabelInterpreter='latex';
axis(FigHandle,[0 TFinal MyMin-MyOffSet MyMax+MyOffSet]);
set(FigHandle,'YTick',floor(MyMin):0.05:ceil(MyMax));
 grid on;



 

 
 if exist('Figures')~=7
    mkdir('Figures'); 
 end

 cd('Figures'); 
  print -dpdf Case39_LQROPF_alphapoint6_LQR_v.pdf
print -depsc2 Case39_LQROPF_alphapoint6_LQR_v
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
plot(t,(ALqrOpf.vVec).','lineWidth',2);
FigHandle=gca; 
set(FigHandle,'box','on');
set(FigHandle,'fontSize',20); 
set(FigHandle,'defaulttextinterpreter','latex');
ytickformat(FigHandle,'%.2f'); 

grid on;
xlabel('Time (sec)', 'FontWeight','bold');
 ylabel('$\mathbf{v}-\mathbf{v}^\mathrm{eq}$ (pu)'); 
  title('Nodal voltages');

FigHandle.TickLabelInterpreter='latex';
axis(FigHandle,[0 TFinal MyMin-MyOffSet MyMax+MyOffSet]);
set(FigHandle,'YTick',floor(MyMin):0.05:ceil(MyMax));
 grid on;



 

 
 if exist('Figures')~=7
    mkdir('Figures'); 
 end

 cd('Figures'); 
  print -dpdf Case39_ALQROPF_alphapoint6_LQR_v.pdf
print -depsc2 Case39_ALQROPF_alphapoint6_LQR_v
cd('..');

%% mech plots


x0=0;
y0=1;
width=8;
height=5;
Figure4=figure('Units','inches',...
'Position',[x0 y0 width height],...
'PaperPositionMode','auto');
set(Figure4, 'Name', 'LQR-OPF:v');
plot(t,(Opf.vVec).','lineWidth',2);
FigHandle=gca; 
set(FigHandle,'box','on');
set(FigHandle,'fontSize',20); 
set(FigHandle,'defaulttextinterpreter','latex');
ytickformat(FigHandle,'%.2f'); 

grid on;
xlabel('Time (sec)', 'FontWeight','bold');
 ylabel('$\mathbf{v}-\mathbf{v}^\mathrm{eq}$ (pu)'); 
  title('Nodal voltages');

FigHandle.TickLabelInterpreter='latex';
axis(FigHandle,[0 TFinal MyMin-MyOffSet MyMax+MyOffSet]);
set(FigHandle,'YTick',floor(MyMin):0.05:ceil(MyMax));
 grid on;


 
 

 
 if exist('Figures')~=7
    mkdir('Figures'); 
 end

 cd('Figures'); 
  print -dpdf Case39_OPF_alphapoint6_LQR_v.pdf
print -depsc2 Case39_OPF_alphapoint6_LQR_v
cd('..');




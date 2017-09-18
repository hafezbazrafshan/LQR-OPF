function  plotsForSolvedCase(casefile )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
% 
% 
% 
% 

load([casefile.DynamicPath, '/',casefile.DynamicFileName]);
cd(DynamicPath);
if exist('figures')~=7
    mkdir('figures');
end
cd('figures'); 
   
get(0, 'DefaultFigureVisible');
set(0, 'DefaultFigureVisible', 'on');
% plots
figx0=0;
figy0=1;
width=8;
height=5;
% 
% 
disp('Plots to be generated');

FreqyMin=min(min(omegaVec))./(2*pi);
FreqyMax=max(max(omegaVec))./(2*pi);
FreqyOffSet=0;
Figure1=figure('Units','inches',...
'Position',[figx0 figy0 width height],...
'PaperPositionMode','auto');
set(Figure1, 'Name', 'GenFreq');
plot(t,omegaVec./(2*pi),'lineWidth',2);
 xlabel('Time (sec)', 'FontWeight','bold');
 ylabel('$\frac{1}{2\pi}\boldmath{\omega}$ (Hz)'); 
axis([0 TFinal FreqyMin-FreqyOffSet FreqyMax+FreqyOffSet]);
set(gca,'box','on');
set(gca,'fontSize',22); 
set(0,'defaulttextinterpreter','latex')
 grid on;
title(' Generator frequencies'); 
  print -dpdf freq.pdf
print -depsc2 freq
% 
% 
% 
% 
% 
% 
AngleyMin=min(min(deltaVec-repmat(deltaS,1,length(t))));
AngleyMax=max(max(deltaVec-repmat(deltaS,1,length(t))));
AngleyOffSet=0;
Figure2=figure('Units','inches',...
'Position',[figx0 figy0 width height],...
'PaperPositionMode','auto');
set(Figure2, 'Name', 'GenAngle');
plot(t,deltaVec-repmat(deltaS,1,length(t)),'lineWidth',2);
 xlabel('Time (sec)', 'FontWeight','bold');
 ylabel('$\boldmath{\delta}-\boldmath{\delta}^s$ (Rad)'); 
axis([0 TFinal AngleyMin-AngleyOffSet AngleyMax+AngleyOffSet]);
set(gca,'box','on');
set(gca,'fontSize',22); 
set(0,'defaulttextinterpreter','latex')
 grid on;
title('Generator angle dev.'); 
  print -dpdf angles.pdf
print -depsc2 angles




eyMin=min(min(eVec-repmat(eS,1,length(t))));
eyMax=max(max(eVec-repmat(eS,1,length(t))));
eyOffSet=0;
Figure3=figure('Units','inches',...
'Position',[figx0 figy0 width height],...
'PaperPositionMode','auto');
set(Figure3, 'Name', 'Gene');
plot(t,eVec-repmat(eS,1,length(t)),'lineWidth',2);
 xlabel('Time (sec)', 'FontWeight','bold');
 ylabel('$\boldmath{e}-\boldmath{e}^s$ (pu)'); 
axis([0 TFinal eyMin-eyOffSet eyMax+eyOffSet]);
set(gca,'box','on');
set(gca,'fontSize',22); 
set(0,'defaulttextinterpreter','latex')
 grid on;
title('Generator EMF dev.'); 
  print -dpdf e.pdf
print -depsc2 e


myMin=min(min(mVec-repmat(mS,1,length(t))));
myMax=max(max(mVec-repmat(mS,1,length(t))));
myOffSet=0;
Figure4=figure('Units','inches',...
'Position',[figx0 figy0 width height],...
'PaperPositionMode','auto');
set(Figure4, 'Name', 'gene');
plot(t,mVec-repmat(mS,1,length(t)),'lineWidth',2);
 xlabel('Time (sec)', 'FontWeight','bold');
 ylabel('$\boldmath{m}-\boldmath{m}^s$ (pu)'); 
axis([0 TFinal myMin-myOffSet myMax+myOffSet]);
set(gca,'box','on');
set(gca,'fontSize',22); 
set(0,'defaulttextinterpreter','latex')
 grid on;
title('Mechanical input dev.'); 
  print -dpdf m.pdf
print -depsc2 m

% 
% 
% 
% 
%     
vyMin=min(min(vVec));
vyMax=max(max(vVec));
vyOffSet=0;
Figure5=figure('Units','inches',...
'Position',[figx0 figy0 width height],...
'PaperPositionMode','auto');
set(Figure5, 'Name', 'Voltage mags');
plot(t,vVec,'lineWidth',2);
 xlabel('Time (sec)', 'FontWeight','bold');
 ylabel('$\boldmath{v}$ (pu)'); 
axis([0 TFinal vyMin-vyOffSet vyMax+vyOffSet]);
set(gca,'box','on');
set(gca,'fontSize',22); 
set(0,'defaulttextinterpreter','latex')
 grid on;
title('Nodal voltage magnitudes'); 
  print -dpdf VoltageMags.pdf
print -depsc2 VoltageMags



thetayMin=min(min(thetaVec));
thetayMax=max(max(thetaVec));
thetayOffSet=0;
Figure6=figure('Units','inches',...
'Position',[figx0 figy0 width height],...
'PaperPositionMode','auto');
set(Figure6, 'Name', 'voltage angles');
plot(t,thetaVec,'lineWidth',2);
 xlabel('Time (sec)', 'FontWeight','bold');
 ylabel('$\boldmath{\theta}$ (Rad)'); 
axis([0 TFinal thetayMin-thetayOffSet thetayMax+thetayOffSet]);
set(gca,'box','on');
set(gca,'fontSize',22); 
set(0,'defaulttextinterpreter','latex')
 grid on;
title('Nodal voltage angles'); 
  print -dpdf VoltageAngles.pdf
print -depsc2 VoltageAngles





pgyMin=min(min(pgVec-repmat(pgS,1,length(t))));
pgyMax=max(max(pgVec-repmat(pgS,1,length(t))));
pgyOffSet=0;
Figure7=figure('Units','inches',...
'Position',[figx0 figy0 width height],...
'PaperPositionMode','auto');
set(Figure7, 'Name', 'voltage angles');
plot(t,pgVec-repmat(pgS,1,length(t)),'lineWidth',2);
 xlabel('Time (sec)', 'FontWeight','bold');
 ylabel('$\boldmath{p}_g-\boldmath{p}_g^s$ (pu)'); 
axis([0 TFinal pgyMin-pgyOffSet pgyMax+pgyOffSet]);
set(gca,'box','on');
set(gca,'fontSize',22); 
set(0,'defaulttextinterpreter','latex')
 grid on;
title('Real  power output dev.'); 
  print -dpdf pg.pdf
print -depsc2 pg



qgyMin=min(min(qgVec-repmat(qgS,1,length(t))));
qgyMax=max(max(qgVec-repmat(qgS,1,length(t))));
qgyOffSet=0;
Figure8=figure('Units','inches',...
'Position',[figx0 figy0 width height],...
'PaperPositionMode','auto');
set(Figure8, 'Name', 'pg');
plot(t,qgVec-repmat(qgS,1,length(t)),'lineWidth',2);
 xlabel('Time (sec)', 'FontWeight','bold');
 ylabel('$\boldmath{q}_g-\boldmath{q}_g^s$ (pu)'); 
axis([0 TFinal qgyMin-qgyOffSet qgyMax+qgyOffSet]);
set(gca,'box','on');
set(gca,'fontSize',22); 
set(0,'defaulttextinterpreter','latex')
 grid on;
title('Reactive power output dev'); 
  print -dpdf qg.pdf
print -depsc2 qg






prefyMin=min(min(prefVec-repmat(prefS,1,length(t))));
prefyMax=max(max(prefVec-repmat(prefS,1,length(t))));
prefyOffSet=0;
Figure9=figure('Units','inches',...
'Position',[figx0 figy0 width height],...
'PaperPositionMode','auto');
set(Figure9, 'Name', 'pg');
plot(t,prefVec-repmat(prefS,1,length(t)),'lineWidth',2);
 xlabel('Time (sec)', 'FontWeight','bold');
 ylabel('$\boldmath{r}- \boldmath{r}^s$ (pu)'); 
axis([0 TFinal prefyMin-prefyOffSet prefyMax+prefyOffSet]);
set(gca,'box','on');
set(gca,'fontSize',22); 
set(0,'defaulttextinterpreter','latex')
 grid on;
title('Governor signal dev.'); 
  print -dpdf pref.pdf
print -depsc2 pref



fyMin=min(min(fVec-repmat(fS,1,length(t))));
fyMax=max(max(fVec-repmat(fS,1,length(t))));
fyOffSet=0;
Figure10=figure('Units','inches',...
'Position',[figx0 figy0 width height],...
'PaperPositionMode','auto');
set(Figure10, 'Name', 'pg');
plot(t,fVec-repmat(fS,1,length(t)),'lineWidth',2);
 xlabel('Time (sec)', 'FontWeight','bold');
 ylabel('$\boldmath{f}-\boldmath{f}^s$ (pu)'); 
axis([0 TFinal fyMin-fyOffSet fyMax+fyOffSet]);
set(gca,'box','on');
set(gca,'fontSize',22); 
set(0,'defaulttextinterpreter','latex')
 grid on;
title('Exciter field voltage dev.'); 
  print -dpdf f.pdf
print -depsc2 f




cd(CurrentDirectory);
close all;
end





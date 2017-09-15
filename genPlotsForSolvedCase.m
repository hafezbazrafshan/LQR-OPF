function  genPlotsForSolvedCase(casefile )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
% 
% 
% 
% 

v2struct(casefile);

if exist('figures')~=7
    mkdir('figures');
end
cd('figures'); 
% plots
figx0=0;
figy0=1;
width=8;
height=5;
% 
% 


freqyMin=min(min(omegaVec))./(2*pi);
freqyMax=max(max(omegaVec))./(2*pi);
freqyOffSet=0;
figure1=figure('Units','inches',...
'Position',[figx0 figy0 width height],...
'PaperPositionMode','auto');
set(figure1, 'Name', 'genfreq');
plot(t,omegaVec./(2*pi),'lineWidth',2);
 xlabel('Time (sec)', 'FontWeight','bold');
 ylabel('$\frac{1}{2\pi}\boldmath{\omega}$ (Hz)'); 
axis([0 Tfinal freqyMin-freqyOffSet freqyMax+freqyOffSet]);
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
angleyMin=min(min(deltaVec));
angleyMax=max(max(deltaVec));
angleyOffSet=0;
figure1=figure('Units','inches',...
'Position',[figx0 figy0 width height],...
'PaperPositionMode','auto');
set(figure1, 'Name', 'genangle');
plot(t,deltaVec,'lineWidth',2);
 xlabel('Time (sec)', 'FontWeight','bold');
 ylabel('$\boldmath{\delta}$ (Rad)'); 
axis([0 Tfinal angleyMin-angleyOffSet angleyMax+angleyOffSet]);
set(gca,'box','on');
set(gca,'fontSize',22); 
set(0,'defaulttextinterpreter','latex')
 grid on;
title('Generator angles'); 
  print -dpdf angles.pdf
print -depsc2 angles




eyMin=min(min(eVec));
eyMax=max(max(eVec));
eyOffSet=0;
figure1=figure('Units','inches',...
'Position',[figx0 figy0 width height],...
'PaperPositionMode','auto');
set(figure1, 'Name', 'gene');
plot(t,eVec,'lineWidth',2);
 xlabel('Time (sec)', 'FontWeight','bold');
 ylabel('$\boldmath{e}$ (pu)'); 
axis([0 Tfinal eyMin-eyOffSet eyMax+eyOffSet]);
set(gca,'box','on');
set(gca,'fontSize',22); 
set(0,'defaulttextinterpreter','latex')
 grid on;
title('Generator EMF'); 
  print -dpdf e.pdf
print -depsc2 e


myMin=min(min(mVec));
myMax=max(max(mVec));
myOffSet=0;
figure1=figure('Units','inches',...
'Position',[figx0 figy0 width height],...
'PaperPositionMode','auto');
set(figure1, 'Name', 'gene');
plot(t,mVec,'lineWidth',2);
 xlabel('Time (sec)', 'FontWeight','bold');
 ylabel('$\boldmath{m}$ (pu)'); 
axis([0 Tfinal myMin-myOffSet myMax+myOffSet]);
set(gca,'box','on');
set(gca,'fontSize',22); 
set(0,'defaulttextinterpreter','latex')
 grid on;
title('Generator mechanical input'); 
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
figure1=figure('Units','inches',...
'Position',[figx0 figy0 width height],...
'PaperPositionMode','auto');
set(figure1, 'Name', 'Voltage mags');
plot(t,vVec,'lineWidth',2);
 xlabel('Time (sec)', 'FontWeight','bold');
 ylabel('$\boldmath{v}$ (pu)'); 
axis([0 Tfinal vyMin-vyOffSet vyMax+vyOffSet]);
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
figure1=figure('Units','inches',...
'Position',[figx0 figy0 width height],...
'PaperPositionMode','auto');
set(figure1, 'Name', 'voltage angles');
plot(t,thetaVec,'lineWidth',2);
 xlabel('Time (sec)', 'FontWeight','bold');
 ylabel('$\boldmath{\theta}$ (Rad)'); 
axis([0 Tfinal thetayMin-thetayOffSet thetayMax+thetayOffSet]);
set(gca,'box','on');
set(gca,'fontSize',22); 
set(0,'defaulttextinterpreter','latex')
 grid on;
title('Nodal voltage angles'); 
  print -dpdf VoltageAngles.pdf
print -depsc2 VoltageAngles





pgyMin=min(min(pgVec));
pgyMax=max(max(pgVec));
pgyOffSet=0;
figure1=figure('Units','inches',...
'Position',[figx0 figy0 width height],...
'PaperPositionMode','auto');
set(figure1, 'Name', 'voltage angles');
plot(t,pgVec,'lineWidth',2);
 xlabel('Time (sec)', 'FontWeight','bold');
 ylabel('$\boldmath{pg}$ (pu)'); 
axis([0 Tfinal pgyMin-pgyOffSet pgyMax+pgyOffSet]);
set(gca,'box','on');
set(gca,'fontSize',22); 
set(0,'defaulttextinterpreter','latex')
 grid on;
title('Generators real electrical power output'); 
  print -dpdf pg.pdf
print -depsc2 pg



qgyMin=min(min(qgVec));
qgyMax=max(max(qgVec));
qgyOffSet=0;
figure1=figure('Units','inches',...
'Position',[figx0 figy0 width height],...
'PaperPositionMode','auto');
set(figure1, 'Name', 'pg');
plot(t,qgVec,'lineWidth',2);
 xlabel('Time (sec)', 'FontWeight','bold');
 ylabel('$\boldmath{qg}$ (pu)'); 
axis([0 Tfinal qgyMin-qgyOffSet qgyMax+qgyOffSet]);
set(gca,'box','on');
set(gca,'fontSize',22); 
set(0,'defaulttextinterpreter','latex')
 grid on;
title('Generators reactive power output'); 
  print -dpdf qg.pdf
print -depsc2 qg






prefyMin=min(min(prefVec));
prefyMax=max(max(prefVec));
prefyOffSet=0;
figure1=figure('Units','inches',...
'Position',[figx0 figy0 width height],...
'PaperPositionMode','auto');
set(figure1, 'Name', 'pg');
plot(t,prefVec,'lineWidth',2);
 xlabel('Time (sec)', 'FontWeight','bold');
 ylabel('$\boldmath{pref}$ (pu)'); 
axis([0 Tfinal prefyMin-prefyOffSet prefyMax+prefyOffSet]);
set(gca,'box','on');
set(gca,'fontSize',22); 
set(0,'defaulttextinterpreter','latex')
 grid on;
title('Generator  governor signal'); 
  print -dpdf pref.pdf
print -depsc2 pref



fyMin=min(min(fVec));
fyMax=max(max(fVec));
fyOffSet=0;
figure1=figure('Units','inches',...
'Position',[figx0 figy0 width height],...
'PaperPositionMode','auto');
set(figure1, 'Name', 'pg');
plot(t,fVec,'lineWidth',2);
 xlabel('Time (sec)', 'FontWeight','bold');
 ylabel('$\boldmath{f}$ (pu)'); 
axis([0 Tfinal fyMin-fyOffSet fyMax+fyOffSet]);
set(gca,'box','on');
set(gca,'fontSize',22); 
set(0,'defaulttextinterpreter','latex')
 grid on;
title('Exciter field voltage'); 
  print -dpdf f.pdf
print -depsc2 f



   
get(0, 'DefaultFigureVisible');
set(0, 'DefaultFigureVisible', 'on');
cd(CurrentDirectory);
end


end


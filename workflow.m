function out=workFlow(CaseFile,StControl,LfControl,Alpha)
%% Guaranteeing that matpower is on the filepath:
% The code uses MATPOWER, and it is important to make sure matpower is
% enabled.
clc;
close all;
CurrentDirectory=pwd;
cd('..'); 
try
cd('matpower6.0/'); 
MatPowerDirectory=pwd;
cd(CurrentDirectory); 
addpath(MatPowerDirectory); 
disp('MATPOWER was sucessfully added to the path'); 
catch 
disp('ERROR: unable to find MATPOWER')
end


% The following is being issued since CVX uses nargin.
%  Matlab has depcrecated nargin and is using nargchk.
% The warnings are annoying so I turned them off.
warning('off','MATLAB:nargchk:deprecated')

global ControlMode
ControlMode=LfControl;

global SteadyStateMode
SteadyStateMode=StControl;


%% Defining some global variables
% Some of these variables refer to the network. For example, Sbase, Ymat
% (the bus admittance matrix), etc.  Some others are specific indices
% within a vector.  Some are LQR parameters.  
% These variables are declared global to avoid extra function
% arguments. 

% system constants [these do not change]
global OMEGAS Sbase N G L NodeSet GenSet LoadSet YMat GMat BMat Cg...
    YffVec YftVec  YtfVec YttVec

%  indices [these  do not change]
global deltaIdx omegaIdx eIdx mIdx  ...
    thetaIdx vIdx pgIdx qgIdx  fIdx prefIdx  SlackIdx GenSlackIdx NonSlackSet GenNonSlackSet

% machine [these do not change]
global  TauVec XdVec XqVec XprimeVec DVec MVec TchVec FreqRVec...
    

% dynamical simulations 
global TFinal TPert FSample NSamples NPertSamples Mass...
PertSet PPertValues QPertValues NoiseVarianceSet 


% initial conditions
global x0 omega0 delta0 e0 m0...
    a0 v0 theta0 pg0 qg0...
    u0 pref0 f0...
    pd0 qd0...
    vg0 thetag0...
    z0 Network


% 0plus conditions
global x0Plus omega0Plus delta0Plus e0Plus m0Plus...
    a0Plus v0Plus theta0Plus pg0Plus qg0Plus...
    u0Plus pref0Plus f0Plus...
    pd0Plus qd0Plus...
    vg0Plus thetag0Plus...
    z0Plus...
    deltaDot0Plus omegaDot0Plus eDot0Plus mDot0Plus


% next time-slot conditions
global xS omegaS deltaS eS mS...
    aS vS thetaS pgS qgS...
    uS prefS fS...
    pdS qdS...
    vgS thetagS...
    zS NetworkS
 
 % global LQR
global  KLQRstep Tlqr


if strcmp(ControlMode,'AGC')
global y0 y0Plus yS yDot0Plus yDotS yIdx...
    ParticipationFactors NumberOfAreas AreaSet TieLineFromSet TieLineToSet...
    ACE0Plus PScheduledS GensPerArea BusesPerArea...
    KI KACE KPG KPflow KSumPG KThetaSlack

end


MatPowerOptions=mpoption('out.all',0); % this suppresses MATPOWER print output
MatPowerOptions = mpoption('model', 'AC', 'pf.alg', 'NR', 'verbose', 0, 'out.all',0); 

%% 1.  Importing the test case
CaseStr=['casefiles/',CaseFile];
Network=loadcase(CaseStr);

% Derive bus admittance matrix and relevant network information
% N is the number of buses
% G is the number of generators
% L is the number of non-generator buses 
% YMat is the complex bus admittance matrix
% GMat is the bus conductance matrix
% BMat is the bus susceptance matrix
% NodeSet is the set of nodes 1:N
% GenSet is the set of generators 
% LoadSet is the set of loads 
% Cg is the generator connection matrix
% YffVec, YftVec, YtfVec, YttVec are branch admittances
[ N,G,L,YMat, GMat, BMat,...
    NodeSet, GenSet, LoadSet,...
   Cg, YffVec, YftVec, YtfVec, YttVec] = networkParams( Network );


% Base values
% (voltage base is not explicitly needed)
% (Power base is required because MatPower stores network powers in 
% actual MVA)
Sbase=Network.baseMVA;
OMEGAS=2*pi*60; % synchronous frequency
SlackIdx=find(Network.bus(:,2)==3);
GenSlackIdx=find(Network.gen(:,1)==SlackIdx);
NonSlackSet=setdiff(LoadSet,SlackIdx);
GenNonSlackSet=setdiff([1:G], GenSlackIdx);

%% Areas:
if strcmp(ControlMode,'AGC')
AreaSet=unique(Network.bus(:,7)); 
NumberOfAreas=length(AreaSet); 
BusesPerArea=cell(NumberOfAreas,1); 
GensPerArea=cell(NumberOfAreas,1);
TieLineFromSet=cell(NumberOfAreas,1);
TieLineToSet=cell(NumberOfAreas,1);
for ii=1:NumberOfAreas
    BusesPerArea{ii,1}=find(Network.bus(:,7)==AreaSet(ii)); 
    GensPerArea{ii,1}=find(Network.bus(gen_set,7)==AreaSet(ii));
    TieLineFromSet{ii,1}= find(and(Network.bus(Network.branch(:,1),7)==AreaSet(ii),Network.bus(Network.branch(:,2),7)~=AreaSet(ii)));
    TieLineToSet{ii,1}=find(and(Network.bus(Network.branch(:,1),7)~=AreaSet(ii),Network.bus(Network.branch(:,2),7)==AreaSet(ii)));
end
end



%% 2.  Populating the initialized steady-state variables from the test case file:
% we first run an opf to ensure problem is feasible with specified limits
% Network.bus(:,3)=1.05*Network.bus(:,3);
% Network.bus(:,4)=1.05*Network.bus(:,4);
% Network.gen(:,9)=1.05*Network.gen(:,9);
% Network.gen(:,4)=1.05*Network.gen(:,4);
Network=runpf(Network,MatPowerOptions);
v0=Network.bus(:,8); % eighth column of bus matrix is voltage magnitude solution
theta0=degrees2radians(Network.bus(:,9)); % nineth column of bus matrix is voltage phase solution
pg0=Network.gen(:,2)./Sbase; % second column of gen matrix is real power set points (or solution for slack bus)
qg0=Network.gen(:,3)./Sbase; % third column of gen matrix is reactive power solutions
pd0=Network.bus(:,3)./Sbase; % third column of bus matrix is real power demand
qd0=Network.bus(:,4)./Sbase; % fourth column of bus matrix is reactive power demand
a0=[v0;theta0;pg0;qg0];

% we increase the maximum limit on active and reactive power generation by
% a little bit to ensure feasibility
Network.gen(:,9)=(max(pg0)+0.00001)*Sbase;
Network.gen(:,4)=(max(qg0)+0.00001)*Sbase;

% Verifying the initial power flow solution:
 [checkpf, checkEqs,realGen_check, reactiveGen_check, ...
    realLoad_check,reactiveLoad_check]=...
   checkPowerFlows(v0,theta0,pg0,qg0, pd0,qd0);
if checkpf==1
    disp('Initial power flow solution was correct'); 
else 
    disp('Initial power flow solution was incorrect'); 
    disp('Check case file'); 
    pause;
end

%% 3.  Adding transient parameters to the network:
% The original MATPOWER case file does not include transient parameters.
% The imported case file has machine constants embedded as mac_con or we add it here.
% We use `mac_con' as a semblance of PST.  
% retrieve machine constants:





if isfield(Network,'mac_con')
Sbase2=Network.mac_con(:,3); 
TauVec=Network.mac_con(:,9);
XdVec=Network.mac_con(:,6).*Sbase./Sbase2;
XqVec=Network.mac_con(:,11).*Sbase./Sbase2;
XprimeVec=Network.mac_con(:,7).*Sbase./Sbase2;
DVec=Network.mac_con(:,17).*Sbase./Sbase2;
MVec=Network.mac_con(:,16)/(pi*60).*Sbase2./Sbase;
% this implementation requires tau_vec, xprime_vec, and xq_vec to be
% nonzero.
TauVec(TauVec==0)=5;
XdVec(XdVec==0)=mean(XdVec(XdVec~=0));
XqVec(XqVec==0)=mean(XqVec(XqVec~=0));
XprimeVec(XprimeVec==0)=mean(XprimeVec(XprimeVec~=0));
clear Sbase2;

else
    
TauVec=repmat(5,G,1);
XdVec=repmat(0.7,G,1);
XqVec=repmat(0.7,G,1);
XprimeVec=repmat(0.1,G,1);
DVec=zeros(G,1);
MVec=0.3*repmat(1,G,1);
    
end

TchVec=0.1*ones(G,1); 
FreqRVec=0.01*ones(G,1).*(2*pi); 


%% 4. Obtain generator internal angles and electromotive force from the power flow solution
% set starting frequency to nominal value
vg0=v0(GenSet);
thetag0=theta0(GenSet);

[ delta0, e0]=obtainGenStates(vg0, thetag0, pg0, qg0 );
omega0=repmat(OMEGAS,G,1); % creating a vector of OMEGAS of size(G,1), for all generator nodes.


%% 5.  Obtaining generator steady-state controls from the power flow solutions and steady-state of states
[m0,f0]=obtainGenControls(delta0,omega0,e0,vg0,thetag0,pg0,qg0, OMEGAS);
x0=[delta0;omega0;e0;m0];
pref0=m0;
if strcmp(ControlMode,'AGC')
y0=zeros(G,1); 
end
u0=[pref0;f0];
%% 6.  Defining the indices of vector z for dynamical simulation:
% *****states*****x
% delta size(G,1)
% omega size(G,1)
% e size(G,1)
deltaIdx=(1:G).';
omegaIdx=(deltaIdx(end)+1:deltaIdx(end)+G).';
eIdx=(omegaIdx(end)+1:omegaIdx(end)+G).';
mIdx=(eIdx(end)+1:eIdx(end)+G).';



%*****algebraic variables*****a
% theta size(N,1)
% v size(N,1)
% pg size(G,1)
% qg size(G,1)
vIdx=(1:N).';
thetaIdx=(vIdx(end)+1:vIdx(end)+N).';
pgIdx=(thetaIdx(end)+1:thetaIdx(end)+G).';
qgIdx=(pgIdx(end)+1:pgIdx(end)+G).';


if strcmp(ControlMode,'AGC')
yIdx=(1:G).';
end

%*****control variables*******u
%pref and f
prefIdx=(1:G).';
fIdx=(prefIdx(end)+1:prefIdx(end)+G).';







%% 7.  Introducing new load for the next OPF time-slot:
disp('Assigning perturbations to load'); 
PRatio=0.1; 
QRatio=0.0484;
PertSet=find(or(Network.bus(:,3)>0, Network.bus(:,4)>0));
PPertValues=PRatio*pd0(PertSet);
QPertValues=QRatio*qd0(PertSet); 
NoiseVarianceSet=0*Network.bus(:,3)/Sbase;
for jj=1:length(PertSet)
   MessageSTR=['Real and reactive power consumption of bus No. ', num2str(PertSet(jj)), ' respectively perturbed by ',...
       num2str(PPertValues(jj)), ' puWatts.', ' and ',  num2str(QPertValues(jj)), ' puVars.']; 
   disp(MessageSTR);
    pause(0.03); 
end

% new steady-state conditions
[pdS,qdS]=loadPert('Steady-State',[],pd0,qd0,PertSet,PPertValues,QPertValues,[],[],[]);

NetworkS=Network;
NetworkS.bus(:,3)= pdS.*Sbase; 
NetworkS.bus(:,4)=qdS.*Sbase; 
 NetworkS.branch(:,[6 7 8])=0; % flow limits are set to zero
 %% 8. Setting LQR parameters
% alpha=0.8;
Tlqr=1000;


%% 9.  Solving the augmentedOPF for the next time-slot
% setting matpower options need in subsequent load-flow

if strcmp(SteadyStateMode,'LQR-OPF')
    % to obtain Vg, pgS, qgS
dpgS=pdS(GenSet)-pd0(GenSet);
dqdgS=qdS(GenSet)-qd0(GenSet);
dpdlS=pdS(LoadSet)-pd0(LoadSet);
dqdlS=qdS(LoadSet)-qd0(LoadSet);
[vgS,pgSNonSlack, thetaSSlack,CompTime,~,~, ~] = ...
    augmentedOPF( delta0, omega0, e0, m0, v0, theta0, pg0, qg0, pref0, f0,...
    NetworkS,...
    dpgS,dqdgS,dpdlS,dqdlS,Alpha); 
NetworkS.gen(:,6)=vgS;
NetworkS.gen(GenNonSlackSet,2)=pgSNonSlack.*Sbase;
NetworkS.bus(NetworkS.bus(:,2)==3,9)=radians2degrees(thetaSSlack);
% run matpower power flow:
[NetworkS,SuccessFlag]=runpf(NetworkS,MatPowerOptions);
elseif strcmp(SteadyStateMode,'OPF')
    TStart=tic;
  [NetworkS, SuccessFlag]=  runopf(NetworkS,MatPowerOptions);
  CompTime=toc(TStart);
end    

%% 10. Obtain a true equilibrium for the next time-slot
vS= NetworkS.bus(:,8);
vgS=vS(GenSet);
thetaS= degrees2radians(NetworkS.bus(:,9));
thetagS=thetaS(GenSet); 
pgS=NetworkS.gen(:,2)./Sbase; 
qgS=NetworkS.gen(:,3)./Sbase;
aS=[vS;thetaS;pgS;qgS];
[ SsCost ] = steadyStateCost(pgS, NetworkS);
    


% check new power flow
 [checkpf2, checkEqs2,realGen_check2, reactiveGen_check2, ...
    realLoad_check2,reactiveLoad_check2]=...
   checkPowerFlows(vS,thetaS,pgS,qgS, pdS,qdS);

if checkpf2==1
    disp('The power flow solution for the second time slot is correct'); 
else 
    disp('The power flow solution for the second time slot in incorrect'); 
    disp('Check the new network conditions and MATPOWER runpf successflag'); 
    pause;
end

[ deltaS,eS]=obtainGenStates(vgS, thetagS, pgS, qgS );
omegaS=repmat(OMEGAS,G,1);

% Obtaining generator steady-state controls
[mS,fS]=obtainGenControls( deltaS,omegaS,eS,vgS,thetagS,pgS,qgS, OMEGAS);
prefS=mS;
xS=[deltaS;omegaS;eS;mS];
uS=[prefS;fS];



%% Run new LQR with setpoints determined to reach the true equilibrium

% [KLQRstep,trCost2,gammaEstimate] = LQRstep( zS, z0,alpha,networkS); the
% output has been verified to be identical to LQRstepCARE (for case9 and
% case39)
% Also, trCost2=gammaEstimate 

[KLQRstep,TrCostEstimate,Gamma] = LQRstepCARE( deltaS, omegaS, eS, mS, ...
    vS, thetaS, pgS, qgS, ...
    prefS, fS,...
    delta0, omega0, e0, m0, ...
    v0, theta0, pg0, qg0, ...
    pref0, f0, ...
    Alpha,NetworkS);

%  12. Define the MASS matrix (The E matrix in $E\dot{x}$ descriptor systems)
if strcmp(ControlMode,'LQR')
    Mass=zeros(length(x0)+length(a0), length(x0)+length(a0));
Mass(sub2ind(size(Mass), [deltaIdx,omegaIdx,eIdx,mIdx],[deltaIdx,omegaIdx,eIdx,mIdx]))=1;
elseif strcmp(ControlMode,'AGC')
Mass=zeros(length(x0)+length(a0)+length(y0), length(x0)+ length(a0)+length(y0));
Mass(sub2ind(size(Mass), [deltaIdx,omegaIdx,eIdx,mIdx],[deltaIdx,omegaIdx,eIdx,mIdx]))=1;
Mass(sub2ind(size(Mass), [length(x0)+length(a0)+yIdx],[length(x0)+length(a0)+yIdx]))=1;
end
 %% 13. Set dynamical simulation parameters:
DynamicSolverOptions = odeset('Mass',Mass,'MassSingular', 'yes','MStateDependence','none', ...
    'RelTol',1e-7,'AbsTol',1e-6,'Stats','off');

TFinal=30;
TPert=0; % NEEDS TO BE CLOSE TO t=0 for increased accuracy 
FSample = 100;
NSamples = TFinal * FSample+1;
NPertSamples = max(TPert,0) * FSample+1;
t = 0:1/FSample:TFinal;
NoiseVector=repmat(NoiseVarianceSet,1,NSamples).*randn(length(NoiseVarianceSet),NSamples);

 
 %% Simulation intial conditions at zero plus
display('Configuring 0plus intial conditions due to disturbance');
[pd0Plus,qd0Plus]=loadPert('Transient',0,pd0,qd0,PertSet, PPertValues,QPertValues, TPert,TFinal,NoiseVector);
delta0Plus=delta0;
omega0Plus=omega0;
e0Plus=e0;
m0Plus=m0;
h1Idx=(1:G).';
h2Idx=(G+1:2*G).';
h3Idx=(2*G+1:3*G).';
h4Idx=(3*G+1:4*G).';
h5Idx=(4*G+1:4*G+L).';
h6Idx=(4*G+L+1:4*G+2*L).';
d0Plus=zeros(h6Idx(end),1);
pdg0Plus=pd0Plus(GenSet);
qdg0Plus=qd0Plus(GenSet);
pdl0Plus=pd0Plus(LoadSet);
qdl0Plus=qd0Plus(LoadSet);
d0Plus(h3Idx)=-pdg0Plus;
d0Plus(h4Idx)=-qdg0Plus;
d0Plus(h5Idx)=-pdl0Plus;
d0Plus(h6Idx)=-qdl0Plus;






if strcmp(ControlMode,'AGC')
     ParticipationFactors=cell(NumberOfAreas,1);
    y0plus=zeros(G,1);
     PscheduledS=zeros(NumberOfAreas,1);
     KI=1;
     KACE=	1;
     KPG=0;
     KSumPG=1;
     KPflow=1;
     KThetaSlack=20;
     for ii=1:NumberOfAreas
          [PFromS,~]=determineLineFlows(TieLineFromSet{ii,1},vS,thetaS);
               [~,PToS]=determineLineFlows(TieLineToSet{ii,1},vS,thetaS); 
               PScheduledS(ii,1)=sum(PFromS)-sum(PToS);
ParticipationFactors{ii,1}=pgS(GensPerArea{ii,1})./sum(pgS(GensPerArea{ii,1}));

     end

[~, ACE0Plus, PMeasured0plus, OmegaMeasured0plus] = agcParams(omega0Plus,y0Plus, v0Plus,theta0Plus, pg0Plus);
end

% % from this point on we would like to have a closed loop system using
% % KLQRstep and the optimal controller set-points 
% % 
% % First
if N<250
InitialOptions= optimoptions('fsolve','Display','Iter','Algorithm','levenberg-marquardt','InitDamping',0.5, 'ScaleProblem','jacobian',...
    'SpecifyObjectiveGradient',true,'MaxIterations',50,'MaxFunctionEvaluations',50,'OptimalityTolerance',1e-6,...
    'PlotFcn',@optimplotfval);

[a0Plus]=fsolve(@( a) hTildeAlgebraicFunctionVectorized(...
    delta0Plus,omega0Plus,e0Plus,m0Plus, ...
    a(vIdx), a(thetaIdx), a(pgIdx), a(qgIdx),...
  d0Plus), [v0;theta0;pg0;qg0],InitialOptions);


v0Plus=a0Plus(vIdx);
theta0Plus=a0Plus(thetaIdx); 
pg0Plus=a0Plus(pgIdx); 
qg0Plus=a0Plus(qgIdx);


InitialOptions= optimoptions('fsolve','Display','Iter','Algorithm','trust-region','InitDamping',0.1, 'ScaleProblem','jacobian',...
    'SpecifyObjectiveGradient',true,'MaxIterations',100000,'MaxFunctionEvaluations',50000,'OptimalityTolerance',1e-6,...
    'PlotFcn',@optimplotfval);
[a0Plus]=fsolve(@( a) hTildeAlgebraicFunctionVectorized(...
    delta0Plus,omega0Plus,e0Plus,m0Plus, ...
    a(vIdx), a(thetaIdx), a(pgIdx), a(qgIdx),...
  d0Plus), a0Plus,InitialOptions);
v0Plus=a0Plus(vIdx);
theta0Plus=a0Plus(thetaIdx); 
pg0Plus=a0Plus(pgIdx); 
qg0Plus=a0Plus(qgIdx);
% 
% 
[ deltaDot0Plus, omegaDot0Plus, eDot0Plus, mDot0Plus ] = gTildeFunctionVectorized(...
    delta0Plus, omega0Plus, e0Plus,m0Plus,...
     v0Plus,theta0Plus,pg0Plus,qg0Plus);
%  
% 
DynamicSolverOptions.Jacobian=@dynamicsJacobian; % this is for the ode solver
Znew0=[delta0Plus; omega0Plus; e0Plus; m0Plus; v0Plus; theta0Plus; pg0Plus; qg0Plus];
ZDot0Plus=zeros(size(Znew0));
ZDot0Plus([deltaIdx;omegaIdx;eIdx;mIdx])=[deltaDot0Plus;omegaDot0Plus;eDot0Plus;mDot0Plus];
DynamicSolverOptions.InitialSlope=ZDot0Plus;
% 
% 
% 
% %%

disp('Running stabilityOPF dynamics'); 
[~,ZNEW]=ode15s(@(t,znew)...
    runDynamics(t,znew,NoiseVector), t, Znew0, DynamicSolverOptions);
ZNEW = transpose(ZNEW);
disp('Finished stabilityOPF dynamics'); 
% 
% 
% 
% 
% 
% 
% 
%  
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% %%
% 
% 
[deltaVec, omegaVec, eVec, mVec,...
    thetaVec, vVec, pgVec, qgVec, ...
    prefVec,fVec, ...
    ploadVec, qloadVec,...
    deltaDotVec, omegaDotVec, eDotVec, mDotVec] =...
    retrieveOutput( t, ZNEW , NoiseVector);
% 
[ SanityCheck1,SanityCheck2,SanityCheck3 , Success] = sanityCheck(...
    deltaVec, omegaVec, eVec, mVec, ...
    thetaVec, vVec, pgVec, qgVec,...
    prefVec, fVec, ...
    ploadVec,qloadVec,...
    deltaDotVec, omegaDotVec, eDotVec, mDotVec);
% 
[ TrCost] = calculateTrCostUsingIntegration(pgS, qgS, Alpha,...
    deltaVec, omegaVec, eVec, mVec, prefVec, fVec, ...
    deltaS, omegaS, eS, mS, prefS, fS);
  
end
% 
% 
if exist('Results')~=7
    mkdir('Results'); 
end

cd('Results'); 
if exist(CaseFile)~=7
    mkdir(CaseFile);
end
% 
cd(CaseFile);

if exist(StControl)~=7
    mkdir(StControl);
end
% 
cd(StControl);
% 
if exist(LfControl)~=7
    mkdir(LfControl);
end
% 
cd(LfControl)
savename=[CaseFile,'_',StControl,'_',LfControl,'alphapoint',num2str(ceil(Alpha*10))];
save(savename); 

outname=['Results/',CaseFile,'/',StControl,'/',LfControl,'/',savename];
out=load(outname);
% 
% 
% 
% 
% 
% if exist('figures')~=7
%     mkdir('figures');
% end
% cd('figures'); 
% % plots
% figx0=0;
% figy0=1;
% width=8;
% height=5;
% % 
% % 
% 
% % set(0,'DefaultFigureVisible','off'); % suppressing pop-up figure windows
% freqyMin=min(min(omegaVec))./(2*pi);
% freqyMax=max(max(omegaVec))./(2*pi);
% freqyOffSet=0;
% figure1=figure('Units','inches',...
% 'Position',[figx0 figy0 width height],...
% 'PaperPositionMode','auto');
% set(figure1, 'Name', 'genfreq');
% plot(t,omegaVec./(2*pi),'lineWidth',2);
%  xlabel('Time (sec)', 'FontWeight','bold');
%  ylabel('$\frac{1}{2\pi}\boldmath{\omega}$ (Hz)'); 
% axis([0 Tfinal freqyMin-freqyOffSet freqyMax+freqyOffSet]);
% set(gca,'box','on');
% set(gca,'fontSize',22); 
% set(0,'defaulttextinterpreter','latex')
%  grid on;
% title(' Generator Frequencies'); 
%   print -dpdf freq.pdf
% print -depsc2 freq
% % 
% % 
% % 
% % 
% % 
% % 
% angleyMin=min(min(deltaVec));
% angleyMax=max(max(deltaVec));
% angleyOffSet=0;
% figure1=figure('Units','inches',...
% 'Position',[figx0 figy0 width height],...
% 'PaperPositionMode','auto');
% set(figure1, 'Name', 'genangle');
% plot(t,deltaVec,'lineWidth',2);
%  xlabel('Time (sec)', 'FontWeight','bold');
%  ylabel('$\boldmath{\delta}$ (Rad)'); 
% axis([0 Tfinal angleyMin-angleyOffSet angleyMax+angleyOffSet]);
% set(gca,'box','on');
% set(gca,'fontSize',22); 
% set(0,'defaulttextinterpreter','latex')
%  grid on;
% title('Generator angles'); 
%   print -dpdf angles.pdf
% print -depsc2 angles
% 
% 
% 
% 
% eyMin=min(min(eVec));
% eyMax=max(max(eVec));
% eyOffSet=0;
% figure1=figure('Units','inches',...
% 'Position',[figx0 figy0 width height],...
% 'PaperPositionMode','auto');
% set(figure1, 'Name', 'gene');
% plot(t,eVec,'lineWidth',2);
%  xlabel('Time (sec)', 'FontWeight','bold');
%  ylabel('$\boldmath{e}$ (pu)'); 
% axis([0 Tfinal eyMin-eyOffSet eyMax+eyOffSet]);
% set(gca,'box','on');
% set(gca,'fontSize',22); 
% set(0,'defaulttextinterpreter','latex')
%  grid on;
% title('Generator angles'); 
%   print -dpdf e.pdf
% print -depsc2 e
% 
% 
% myMin=min(min(mVec));
% myMax=max(max(mVec));
% myOffSet=0;
% figure1=figure('Units','inches',...
% 'Position',[figx0 figy0 width height],...
% 'PaperPositionMode','auto');
% set(figure1, 'Name', 'gene');
% plot(t,mVec,'lineWidth',2);
%  xlabel('Time (sec)', 'FontWeight','bold');
%  ylabel('$\boldmath{m}$ (pu)'); 
% axis([0 Tfinal myMin-myOffSet myMax+myOffSet]);
% set(gca,'box','on');
% set(gca,'fontSize',22); 
% set(0,'defaulttextinterpreter','latex')
%  grid on;
% title('Generator angles'); 
%   print -dpdf m.pdf
% print -depsc2 m
% 
% % 
% % 
% % 
% % 
% %     
% vyMin=min(min(vVec));
% vyMax=max(max(vVec));
% vyOffSet=0;
% figure1=figure('Units','inches',...
% 'Position',[figx0 figy0 width height],...
% 'PaperPositionMode','auto');
% set(figure1, 'Name', 'Voltage mags');
% plot(t,vVec,'lineWidth',2);
%  xlabel('Time (sec)', 'FontWeight','bold');
%  ylabel('$\boldmath{v}$ (pu)'); 
% axis([0 Tfinal vyMin-vyOffSet vyMax+vyOffSet]);
% set(gca,'box','on');
% set(gca,'fontSize',22); 
% set(0,'defaulttextinterpreter','latex')
%  grid on;
% title('Generator angles'); 
%   print -dpdf VoltageMags.pdf
% print -depsc2 VoltageMags
% 
% 
% 
% thetayMin=min(min(thetaVec));
% thetayMax=max(max(thetaVec));
% thetayOffSet=0;
% figure1=figure('Units','inches',...
% 'Position',[figx0 figy0 width height],...
% 'PaperPositionMode','auto');
% set(figure1, 'Name', 'voltage angles');
% plot(t,thetaVec,'lineWidth',2);
%  xlabel('Time (sec)', 'FontWeight','bold');
%  ylabel('$\boldmath{\theta}$ (Rad)'); 
% axis([0 Tfinal thetayMin-thetayOffSet thetayMax+thetayOffSet]);
% set(gca,'box','on');
% set(gca,'fontSize',22); 
% set(0,'defaulttextinterpreter','latex')
%  grid on;
% title('Generator angles'); 
%   print -dpdf VoltageAngles.pdf
% print -depsc2 VoltageAngles
% 
% 
% 
% 
% 
% pgyMin=min(min(pgVec));
% pgyMax=max(max(pgVec));
% pgyOffSet=0;
% figure1=figure('Units','inches',...
% 'Position',[figx0 figy0 width height],...
% 'PaperPositionMode','auto');
% set(figure1, 'Name', 'voltage angles');
% plot(t,pgVec,'lineWidth',2);
%  xlabel('Time (sec)', 'FontWeight','bold');
%  ylabel('$\boldmath{pg}$ (pu)'); 
% axis([0 Tfinal pgyMin-pgyOffSet pgyMax+pgyOffSet]);
% set(gca,'box','on');
% set(gca,'fontSize',22); 
% set(0,'defaulttextinterpreter','latex')
%  grid on;
% title('Generated power'); 
%   print -dpdf pg.pdf
% print -depsc2 pg
% 
% 
% 
% qgyMin=min(min(qgVec));
% qgyMax=max(max(qgVec));
% qgyOffSet=0;
% figure1=figure('Units','inches',...
% 'Position',[figx0 figy0 width height],...
% 'PaperPositionMode','auto');
% set(figure1, 'Name', 'pg');
% plot(t,qgVec,'lineWidth',2);
%  xlabel('Time (sec)', 'FontWeight','bold');
%  ylabel('$\boldmath{qg}$ (pu)'); 
% axis([0 Tfinal qgyMin-qgyOffSet qgyMax+qgyOffSet]);
% set(gca,'box','on');
% set(gca,'fontSize',22); 
% set(0,'defaulttextinterpreter','latex')
%  grid on;
% title('Reactive powers'); 
%   print -dpdf qg.pdf
% print -depsc2 qg
% 
% 
% 
% 
% 
% 
% prefyMin=min(min(prefVec));
% prefyMax=max(max(prefVec));
% prefyOffSet=0;
% figure1=figure('Units','inches',...
% 'Position',[figx0 figy0 width height],...
% 'PaperPositionMode','auto');
% set(figure1, 'Name', 'pg');
% plot(t,prefVec,'lineWidth',2);
%  xlabel('Time (sec)', 'FontWeight','bold');
%  ylabel('$\boldmath{pref}$ (pu)'); 
% axis([0 Tfinal prefyMin-prefyOffSet prefyMax+prefyOffSet]);
% set(gca,'box','on');
% set(gca,'fontSize',22); 
% set(0,'defaulttextinterpreter','latex')
%  grid on;
% title('Generator reference'); 
%   print -dpdf pref.pdf
% print -depsc2 pref
% 
% 
% 
% fyMin=min(min(fVec));
% fyMax=max(max(fVec));
% fyOffSet=0;
% figure1=figure('Units','inches',...
% 'Position',[figx0 figy0 width height],...
% 'PaperPositionMode','auto');
% set(figure1, 'Name', 'pg');
% plot(t,fVec,'lineWidth',2);
%  xlabel('Time (sec)', 'FontWeight','bold');
%  ylabel('$\boldmath{f}$ (pu)'); 
% axis([0 Tfinal fyMin-fyOffSet fyMax+fyOffSet]);
% set(gca,'box','on');
% set(gca,'fontSize',22); 
% set(0,'defaulttextinterpreter','latex')
%  grid on;
% title('Generator emf'); 
%   print -dpdf f.pdf
% print -depsc2 f
% 
% 
% 
%    
% get(0, 'DefaultFigureVisible');
% set(0, 'DefaultFigureVisible', 'on');
cd(CurrentDirectory);
% end
% 


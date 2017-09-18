function [ out ] = steadyStateDriver(CaseFile,SsControl, Alpha, Tlqr, PRatio, QRatio )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here



CurrentDirectory=pwd;
cd('..');
disp('Configuring MATPOWER'); 
PreSuccessStr=['.........................................................................'];
PauseTime=0.5;
try
cd('matpower6.0/'); 
MatPowerDirectory=pwd;
cd(CurrentDirectory); 
addpath(MatPowerDirectory); 
disp([PreSuccessStr,'Successful']); 
pause(PauseTime);
catch 
disp('ERROR: unable to find MATPOWER')
pause;
end




% The following is being issued since CVX uses nargin.
%  Matlab has depcrecated nargin and is using nargchk.
% The warnings are annoying so I turned them off.
warning('off','MATLAB:nargchk:deprecated')





global SteadyStateMode
SteadyStateMode=SsControl;


%% Defining some global variables
% Some of these variables refer to the network. For example, Sbase, Ymat
% (the bus admittance matrix), etc.  Some others are specific indices
% within a vector.  Some are LQR parameters.  
% These variables are declared global to avoid extra function
% arguments. 

% system constants [these do not change]
global OMEGAS Sbase N G L NodeSet GenSet LoadSet NodeLabels GenLabels LoadLabels...
    YMat GMat BMat Cg...
    YffVec YftVec  YtfVec YttVec

%  indices [these  do not change]
global deltaIdx omegaIdx eIdx mIdx  ...
    thetaIdx vIdx pgIdx qgIdx  fIdx prefIdx  SlackIdx GenSlackIdx NonSlackSet GenNonSlackSet

% machine [these do not change]
global  TauVec XdVec XqVec XprimeVec DVec MVec TchVec FreqRVec...
    

% dynamical simulations 
global ...
PertSet PPertValues QPertValues 


% initial conditions
global x0 omega0 delta0 e0 m0...
    a0 v0 theta0 pg0 qg0...
    u0 pref0 f0...
    pd0 qd0...
    vg0 thetag0...
    z0 Network




% next time-slot conditions
global xS omegaS deltaS eS mS...
    aS vS thetaS pgS qgS...
    uS prefS fS...
    pdS qdS...
    vgS thetagS...
    zS NetworkS
 
 % global LQR
global  KLQRstep 




MatPowerOptions=mpoption('out.all',0); % this suppresses MATPOWER print output
MatPowerOptions = mpoption('model', 'AC', 'pf.alg', 'NR', 'verbose', 0, 'out.all',0); 

%% 1.  Importing the test case
disp(['Importing ', CaseFile]);
CaseStr=['casefiles/',CaseFile];
Network=loadcase(CaseStr);
disp([PreSuccessStr,'Successful']);
pause(PauseTime);


disp('Populating steady-state network parameters');

[ N,G,L,YMat, GMat, BMat,...
    NodeSet, GenSet, LoadSet,NodeLabels,GenLabels,LoadLabels,...
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


disp([PreSuccessStr,'Successful']);
pause(PauseTime);



disp('Running initial load flow to obtain initial algebraic variables'); 
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




if sum(pg0.*Sbase>=Network.gen(:,9))>0
    Index=find(pg0.*Sbase>=Network.gen(:,9));
    Network.gen(Index,9)=1.01*pg0(Index)*Sbase;
end

if sum(qg0.*Sbase>=Network.gen(:,4))>0
    Index=find(qg0.*Sbase>Network.gen(:,4));
    Network.gen(Index,4)=1.01*qg0(Index)*Sbase;
end


% Verifying the initial power flow solution:
 [checkpf, checkEqs,realGen_check, reactiveGen_check, ...
    realLoad_check,reactiveLoad_check]=...
   checkPowerFlows(v0,theta0,pg0,qg0, pd0,qd0);
if checkpf==1
   disp([PreSuccessStr,'Successful']);
pause(PauseTime);
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



disp('Populating dynamic machine parameters');


if isfield(Network,'mac_con')
    disp('Machine data available');

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
        disp('Machine data not available, Synthetic data is used');

TauVec=repmat(5,G,1);
XdVec=repmat(0.7,G,1);
XqVec=repmat(0.3,G,1);
XprimeVec=repmat(0.06,G,1);
DVec=zeros(G,1);
MVec=0.3*repmat(1,G,1);
    
end

TchVec=0.2*ones(G,1); 
FreqRVec=0.02*ones(G,1).*(2*pi); 

   disp([PreSuccessStr,'Successful']);
pause(PauseTime);

%% 4. Obtain generator internal angles and electromotive force from the power flow solution
% set starting frequency to nominal value
disp('Determining initial machine states from initial algebraic values');
vg0=v0(GenSet);
thetag0=theta0(GenSet);

[ delta0, e0]=obtainGenStates(vg0, thetag0, pg0, qg0 );
omega0=repmat(OMEGAS,G,1); % creating a vector of OMEGAS of size(G,1), for all generator nodes.
   disp([PreSuccessStr,'Successful']);
pause(PauseTime);


%% 5.  Obtaining generator steady-state controls from the power flow solutions and steady-state of states
disp('Determining intial control inputs from initial load flow and state values');
[m0,f0]=obtainGenControls(delta0,omega0,e0,vg0,thetag0,pg0,qg0, OMEGAS);
x0=[delta0;omega0;e0;m0];
pref0=m0;

u0=[pref0;f0];
   disp([PreSuccessStr,'Successful']);
pause(PauseTime);


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


%*****control variables*******u
%pref and f
prefIdx=(1:G).';
fIdx=(prefIdx(end)+1:prefIdx(end)+G).';


%% 7.  Introducing new load for the next OPF time-slot:
disp('Assigning perturbations to load'); 
PertSet=find(or(Network.bus(:,3)>0, Network.bus(:,4)>0));
PPertValues=PRatio*pd0(PertSet);
QPertValues=QRatio*qd0(PertSet); 
% for jj=1:length(PertSet)
% %    MessageSTR=['Modified (Pd,Qd) Bus ', num2str(PertSet(jj)), ' by ',...
% %        num2str(PPertValues(jj)), '+j',  num2str(QPertValues(jj)), 'pu']; 
%    disp(MessageSTR);
%     pause(0.3); 
% end
disp(['Modifying (Pd,Qd) at all buses by ', num2str(PRatio*100), 'Percent with PF=0.9']);

% new steady-state conditions
[pdS,qdS]=loadPert('Steady-State',[],pd0,qd0,PertSet,PPertValues,QPertValues,[],[],[]);

NetworkS=Network;
NetworkS.bus(:,3)= pdS.*Sbase; 
NetworkS.bus(:,4)=qdS.*Sbase; 
   disp([PreSuccessStr,'Successful']);
pause(PauseTime);
 NetworkS.branch(:,[6 7 8])=0; % flow limits are set to zero

 
 %% 9.  Solving the augmentedOPF for the next time-slot
% setting matpower options need in subsequent load-flow
disp(['Steady state optimization requested is ', SteadyStateMode]);
disp(['Running ', SteadyStateMode]);
switch SteadyStateMode
    case 'OPF'
            TStart=tic;
  [NetworkS, SuccessFlag]=  runopf(NetworkS,MatPowerOptions);
  [NetworkS,SuccessFlag]=runpf(NetworkS,MatPowerOptions);
SsObjEst=[];
  CompTime=toc(TStart);      
  
    case 'LQR-OPF'
               TStart=tic;
         [vgS,pgS, thetaSSlack,SsObjEst, ~] = ...
              LQROPF( delta0, omega0, e0, m0, v0, theta0, pg0, qg0, pref0, f0,...
    NetworkS,...
   pdS,qdS,pd0,qd0,Alpha,Tlqr); 
%% SECTION TITLE
% DESCRIPTIVE TEXT
NetworkS.gen(:,6)=vgS;
NetworkS.gen(GenNonSlackSet,2)=pgS(GenNonSlackSet).*Sbase;
NetworkS.bus(NetworkS.bus(:,2)==3,9)=radians2degrees(thetaSSlack);
[NetworkS,SuccessFlag]=runpf(NetworkS,MatPowerOptions);
         CompTime=toc(TStart);
    case 'ALQR-OPF'
          TStart=tic;
            [vgS,pgS, thetaSSlack,SsObjEst, ~] = ...
              ALQROPF( delta0, omega0, e0, m0, v0, theta0, pg0, qg0, pref0, f0,...
    NetworkS,...
    pdS,qdS,pd0,qd0,Alpha,Tlqr); 
        NetworkS.gen(:,6)=vgS;
NetworkS.gen(GenNonSlackSet,2)=pgS(GenNonSlackSet).*Sbase;
NetworkS.bus(NetworkS.bus(:,2)==3,9)=radians2degrees(thetaSSlack);
[NetworkS,SuccessFlag]=runpf(NetworkS,MatPowerOptions);
  CompTime=toc(TStart);
    case 'DLQR-OPF'


end
% run matpower power flow:
  if SuccessFlag==1
         disp([PreSuccessStr,'Successful']);
               disp(['Steady state optimization ', SteadyStateMode, ' took ', num2str(CompTime), ' Seconds']);
pause(PauseTime);
  else
      disp([PreSuccessStr,'Failed!!']);
               disp(['Steady state optimization ', SteadyStateMode, ' Failed']);
pause;
  end

%% 10. Obtain a true equilibrium for the next time-slot
disp(['Retrieiving new steady-state algebraic variables']); 
vS= NetworkS.bus(:,8);
vgS=vS(GenSet);
thetaS= degrees2radians(NetworkS.bus(:,9));
thetagS=thetaS(GenSet); 
pgS=NetworkS.gen(:,2)./Sbase; 
qgS=NetworkS.gen(:,3)./Sbase;
% aS=[vS;thetaS;pgS;qgS];
[ SsCost ] = steadyStateCost(pgS, NetworkS);
          disp([PreSuccessStr,'Successful']);
   pause(PauseTime);



% check new power flow
    disp('Checking whether new steady-state satisfies load flow'); 
 [checkpf2, checkEqs2,realGen_check2, reactiveGen_check2, ...
    realLoad_check2,reactiveLoad_check2]=...
   checkPowerFlows(vS,thetaS,pgS,qgS, pdS,qdS);

if checkpf2==1

        disp([PreSuccessStr,'Successful']);
   pause(PauseTime);
else 
    disp('The power flow solution for the second time slot in incorrect'); 
    disp('Check the new network conditions and MATPOWER runpf successflag'); 
    pause;
end



disp('Determining new machine states from new algebraic values');
[ deltaS,eS]=obtainGenStates(vgS, thetagS, pgS, qgS );
omegaS=repmat(OMEGAS,G,1);

% Obtaining generator steady-state controls
[mS,fS]=obtainGenControls( deltaS,omegaS,eS,vgS,thetagS,pgS,qgS, OMEGAS);
prefS=mS;



%% Run new LQR with setpoints determined to reach the true equilibrium
% 
% [KLQRstep,TrCostEstimate,Gamma] = LQRstep( deltaS, omegaS, eS, mS, ...
%     vS, thetaS, pgS, qgS, ...
%     prefS, fS,...
%     delta0, omega0, e0, m0, ...
%     v0, theta0, pg0, qg0, ...
%     pref0, f0, ...
%     Alpha,NetworkS);

[KLQRstep,TrCostEstimate,Gamma, Asys, Bsys] = LQRstepCARE( deltaS, omegaS, eS, mS, ...
    vS, thetaS, pgS, qgS, ...
    prefS, fS,...
    delta0, omega0, e0, m0, ...
    v0, theta0, pg0, qg0, ...
    pref0, f0, ...
    Alpha,Tlqr,NetworkS);

TotalCostEstimate=SsCost+TrCostEstimate;

MaxEigen=max(real(eig(full(Asys))));
if MaxEigen<0.00
    disp('System is automatically stable'); 
    pause(PauseTime);
else
    disp(['Max. eigen value is ', num2str(MaxEigen)]);
    pause(PauseTime);
end


%% Saving steady-state optimization
if exist('Results')~=7
    mkdir('Results'); 
end

cd('Results'); 
if exist(CaseFile)~=7
    mkdir(CaseFile);
end
% 
cd(CaseFile);

if exist(SsControl)~=7
    mkdir(SsControl);
end
% 
cd(SsControl);
SteadyStateFileName=[CaseFile,'_',SsControl,'_','alphapoint',num2str(ceil(Alpha*10))];
SteadyStatePath=pwd;
save(SteadyStateFileName); 
out=load(SteadyStateFileName);
cd(CurrentDirectory);


end


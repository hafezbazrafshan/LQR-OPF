function [vgS,pgS, thetaSSlack,SsObjEst, Gamma,K...
  RiccatiSVec,KVec, SsCostVec,GammaVec,TrCostEstimateVec,ObjValue, ItSuccess] = ...
   ALQROPF( delta0, omega0, e0, m0,...
    v0,theta0, pg0, qg0, ...
    pref0, f0,...
    NetworkS,...
    pdS,qdS,pd0,qd0, Alpha)
%AUGMENTEDOPF  implements augmented opf per equation (18) CDC 2016. 
%    [vs,thetas, pgs,qgs, K ] = augmentedOPF( z0,...
 %   deltaploadg,deltaqloadg,deltaploadl,deltaqloadl) implements the
 %   augmetned OPF based on linear approximation of a known equilibrium z0
 % for the power systems described by nonlinear equations g(x,a,u) and
 % h(x,a,u). 
 %
 % Description of Outputs: 
 % 1. vs: the calculated optimal steady-state voltage magnitude, size(N,1).
 % 2. thetas: the calculated optimal steady-state voltage angle in radians,
 % size(N,1). 
 % 3. pgs: the calculated optimal steady-state real power injection
 % (setpoints) in pu Watts, size(G,1). 
 % 4. qgs: the calculated optimal steady-state reactive power injection 
% in pu Vars, size(G,1). 
% 5. K: the calculated optimal linear feedback gain, size(2*G,3*G)
% 6. ssCost: is the calculated steady-state cost of real power generation
% 7. trCostEstimate: is the gama---an estimate of the transient cost
% 
% Description of Inputs: 
% 1. z0: the  equilibrium point used for linearization
% 2. deltaploadg: the difference of the new desired  real power to the initial load level
% for generator nodes, size(G,1).
% 3. deltaqloadg: the difference of the new desired  reactive power to the initial load level
% for generator nodes, size(G,1).
% 4. deltaploadl: the difference of the new desired  real power to the initial load level
% for load nodes, size(L,1).
% 5. deltaqloadl: the difference of the new desired  reactive power to the initial load level
% for generator nodes, size(L,1)
% 6. alpha: alpha
% 7. Tlqr: is the Tlqr factor for transient control
% See also approxOPF, LQRstep
%
% Required:
% 


% system constants [these do not change]
global OMEGAS Sbase N G L NodeSet GenSet LoadSet YMat GMat BMat Cg...
    YffVec YftVec  YtfVec YttVec

%  indices [these  do not change]
global deltaIdx omegaIdx eIdx mIdx  ...
    thetaIdx vIdx pgIdx qgIdx  fIdx prefIdx  SlackIdx GenNonSlackSet

% machine [these do not change]
global  TauVec XdVec XqVec XprimeVec DVec MVec TchVec FreqRVec...


global Tlqr

 
    %% Obtaining the jacobians:

 [ gx,ga,gu ] = gFunctionJacobVectorized(delta0, omega0, e0, m0, v0, theta0, pg0, qg0);
[ hx, ha] = hFunctionJacobVectorized( delta0, omega0, e0, m0, v0, theta0, pg0, qg0);

haInv=mldivide(ha,speye(size(ha)));
Asys=gx-ga*haInv*hx;
Bsys=gu;






%%
x0=[delta0;omega0;e0;m0];
a0=[v0;theta0;pg0;qg0];
u0=[pref0;f0];


%% disturbance
h1Idx=1:G;
h2Idx=G+1:2*G;
h3Idx=2*G+1:3*G;
h4Idx=3*G+1:4*G;
h5Idx=4*G+1:4*G+L;
h6Idx=4*G+L+1:4*G+2*L;
deltaD=zeros(h6Idx(end),1);


dpgS=pdS(GenSet)-pd0(GenSet);
dqdgS=qdS(GenSet)-qd0(GenSet);
dpdlS=pdS(LoadSet)-pd0(LoadSet);
dqdlS=qdS(LoadSet)-qd0(LoadSet);

deltaD(h3Idx)=-dpgS;
deltaD(h4Idx)=-dqdgS;
deltaD(h5Idx)=-dpdlS;
deltaD(h6Idx)=-dqdlS;




%% Alternate project-minimize algorithm
 % setting up maximum number of iterations and initial objectives
MaxIt=2;
ObjValue=-inf+zeros(MaxIt,1); 
ZsVec=sparse(length([x0;a0;u0]),MaxIt);
RiccatiSVec=sparse(16*(G^2),MaxIt);
SsCostVec=sparse(1,MaxIt); 
TrCostEstimateVec=sparse(1,MaxIt);
GammaVec=sparse(1,MaxIt);
KVec=sparse(2*G*4*G,MaxIt);



% solve CARE at z0
[Qinv,Rinv]=QinvRinv(pg0,qg0,Alpha, speye(4*G), speye(2*G),NetworkS); 
[RiccatiS,EigValues,FeedBackGain,Report]=care(Asys,Bsys,mldivide(Qinv,speye(size(Qinv))),mldivide(Rinv,speye(size(Rinv))));



%solve OPF with P0
cvx_begin quiet
cvx_solver SDPT3
variables xs(4*G,1) as(2*N+2*G,1) us(2*G,1) 


deltas=xs(deltaIdx); 
omegas=xs(omegaIdx); 
es=xs(eIdx); 
ms=xs(mIdx);
vs=as(vIdx);
thetas=as(thetaIdx); 
pgs=as(pgIdx); 
qgs=as(qgIdx); 
prefs=us(prefIdx); 
fs=us(fIdx);


SsCost=steadyStateCost(pgs,NetworkS);

minimize( SsCost+ (Tlqr/2)*quad_form((xs-x0),RiccatiS))
subject to:
omegas==omega0;
zeros(4*G,1) == gx*(xs-x0)+ga*(as-a0)+gu*(us-u0);
deltaD==hx*(xs-x0)+ha*(as-a0);
NetworkS.bus(:,13)<= vs<=NetworkS.bus(:,12);
NetworkS.gen(:,5)./Sbase<=qgs<=NetworkS.gen(:,4)./Sbase;
NetworkS.gen(:,10)./Sbase <= pgs <= NetworkS.gen(:,9)./Sbase; 



 % finds the slack bus
thetas(SlackIdx)==theta0(SlackIdx);
% vs(SlackIdx)==v0(SlackIdx);



cvx_end


% solve CARE at zS
[Qinv,Rinv]=QinvRinv(pgs,qgs,Alpha, speye(4*G), speye(2*G),NetworkS); 
[RiccatiS,EigValues,FeedBackGain,Report]=care(Asys,Bsys,mldivide(Qinv,speye(size(Qinv))),mldivide(Rinv,speye(size(Rinv))));




ZsVec(:,1)=[xs;as;us];
RiccatiSVec(:,1)=RiccatiS(:);
KVec(:,1)=-FeedBackGain(:);
SsCostVec(:,1)=steadyStateCost(pgs,NetworkS);
GammaVec(:,1)=quad_form((xs-x0),RiccatiS);
TrCostEstimateVec(:,1)=(Tlqr/2)*GammaVec(:,1);
ObjValue(1)=SsCostVec(:,1)+TrCostEstimateVec(:,1);
ItSuccess=1;
BestIt=1;



  K=-FeedBackGain;
SsCost=steadyStateCost(pgs,NetworkS);
Gamma=quad_form((xs-x0),RiccatiS);
TrCostEstimate=(Tlqr/2)*Gamma;
SsObjEst=SsCost+TrCostEstimate;

vgS=vs(GenSet);
% pgSNonSlack=pgs(GenNonSlackSet);
pgS=pgs;
thetaSSlack=thetas(SlackIdx);




str=['ObjValue of Alternate minimization at Iter No. ', num2str(1), ' is ', num2str(ObjValue(1)),'\n'];
fprintf(str);

if Alpha ~=0

for ItNo=2:MaxIt


cvx_begin quiet
cvx_solver sedumi
variables xs(4*G,1) as(2*N+2*G,1) us(2*G,1) 


deltas=xs(deltaIdx); 
omegas=xs(omegaIdx); 
es=xs(eIdx); 
ms=xs(mIdx);
vs=as(vIdx);
thetas=as(thetaIdx); 
pgs=as(pgIdx); 
qgs=as(qgIdx); 
prefs=us(prefIdx); 
fs=us(fIdx);


SsCost=steadyStateCost(pgs,NetworkS);

minimize( SsCost+ (Tlqr/2)*quad_form((xs-x0),RiccatiS))
subject to:
omegas==omega0;
% -pi<=thetas<=pi;
zeros(4*G,1) == gx*(xs-x0)+ga*(as-a0)+gu*(us-u0);
deltaD==hx*(xs-x0)+ha*(as-a0);
NetworkS.bus(:,13)<= vs<=NetworkS.bus(:,12);
NetworkS.gen(:,5)./Sbase<=qgs<=NetworkS.gen(:,4)./Sbase;
NetworkS.gen(:,10)./Sbase <= pgs <= NetworkS.gen(:,9)./Sbase; 

 % finds the slack bus
thetas(SlackIdx)==theta0(SlackIdx);
% vs(SlackIdx)==v0(SlackIdx);


cvx_end


% solve CARE at previous zs
[Qinv,Rinv]=QinvRinv(pgs,qgs,Alpha, speye(4*G), speye(2*G),NetworkS); 
[RiccatiS,EigValues,FeedBackGain,Report]=care(Asys,Bsys,mldivide(Qinv,speye(size(Qinv))),mldivide(Rinv,speye(size(Rinv))));
%solve OPF with P0

% optional outputs of augmentedOPFAlternateMinimization
ZsVec(:,ItNo)=[xs;as;us];
RiccatiSVec(:,ItNo)=RiccatiS(:);
KVec(:,ItNo)=-FeedBackGain(:);
SsCostVec(:,ItNo)=steadyStateCost(pgs,NetworkS);
GammaVec(:,ItNo)=quad_form((xs-x0),RiccatiS);
TrCostEstimateVec(:,ItNo)=(Tlqr/2)*GammaVec(:,ItNo);
ObjValue(ItNo)=SsCostVec(:,ItNo)+TrCostEstimateVec(:,ItNo);

str=['ObjValue of Alternate minimization at Iter No. ', num2str(ItNo), ' is ', num2str(ObjValue(ItNo)),'\n'];
fprintf(str);

% display progress


if ObjValue(ItNo)<ObjValue(BestIt)
  BestIt=ItNo;
  K=-FeedBackGain;

SsCost=steadyStateCost(pgs,NetworkS);
Gamma=quad_form((xs-x0),RiccatiS);
TrCostEstimate=(Tlqr/2)*Gamma;
SsObjEst=SsCost+TrCostEstimate;

vgS=vs(GenSet);
% pgSNonSlack=pgs(GenNonSlackSet);
pgS=pgs;
thetaSSlack=thetas(SlackIdx);
end




end
end



% figureAlt=figure('Units','inches',...
%  'Position',[0 1 8 4],...
%  'PaperPositionMode','auto');
%  set(figureAlt, 'Name', 'Alt. Min. Progress');
% plot(1:1:MaxIt, ObjValue(1:1:ItSuccess),'o--');
%   xlabel('It.', 'FontWeight','bold');
%   ylabel('Obj. Value'); 
%  
%  set(gca, 'XTick',[1:1:MaxIt]);
%  set(gca,'box','on');
%  set(gca,'fontSize',22); 
%  set(gca,'defaulttextinterpreter','latex')
%  grid on;
%  title('Alt. Min. Progress'); 


end


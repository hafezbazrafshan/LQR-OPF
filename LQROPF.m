function [vgS,pgS, thetaSSlack, SsObjEst,ModelingTime, Gamma,K ] = LQROPF(...
    delta0, omega0, e0, m0,...
    v0,theta0, pg0, qg0, ...
    pref0, f0,...
    NetworkS,...
    pdS,qdS,pd0,qd0, Alpha,Tlqr)
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




    
%% Obtaining the jacobians:
 [ gx,ga,gu ] = gFunctionJacobVectorized(delta0, omega0, e0, m0, v0, theta0, pg0, qg0);
[ hx, ha] = hFunctionJacobVectorized( delta0, omega0, e0, m0, v0, theta0, pg0, qg0);

haInv=mldivide(ha,speye(size(ha)));
Asys=gx-ga*haInv*hx;
Bsys=gu;



x0=[delta0;omega0;e0;m0];
a0=[v0;theta0;pg0;qg0];
u0=[pref0;f0];

%% costs:









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


cvx_tic
cvx_begin quiet
cvx_solver SDPT3
variables xs(4*G,1) as(2*N+2*G,1) us(2*G,1)  Gamma
variable P(4*G,4*G) symmetric
variable Y(2*G,4*G) 
expression Qinv(4*G,4*G)
expression Rinv(2*G,2*G) 


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
minimize( SsCost+ (Tlqr/2)*Gamma) ;
subject to:










[Qinv,Rinv]=QinvRinv(pgs,qgs,Alpha, Qinv, Rinv,NetworkS);

Qinv=sparse(Qinv);
Rinv=sparse(Rinv);


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



-[-Gamma, xs.'-x0.'; xs-x0, -P] ==semidefinite(4*G+1);
% 
 -[Asys*P + P*Asys.' + Bsys*Y + Y.' * Bsys.', P , Y.';
    P, -Qinv, zeros(4*G, 2*G); 
    Y, zeros(2*G, 4*G), -Rinv] == semidefinite(10*G);
    

cvx_end
TimeElapse=cvx_toc;
ModelingTime=TimeElapsed(3);


K=-Rinv*Bsys.'*conj(inv(P));


vgS=vs(GenSet);
% pgSNonSlack=pgs(GenNonSlackSet);
pgS=pgs;
thetaSSlack=thetas(SlackIdx);
TrCostEstimate=(Tlqr/2)*Gamma;
SsObjEst=SsCost+TrCostEstimate;

end


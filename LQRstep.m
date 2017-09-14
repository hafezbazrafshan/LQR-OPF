function [K, TrCostEstimate,Gamma] = LQRstep( deltaS, omegaS, eS, mS, ...
    vS, thetaS, pgS, qgS, ...
    prefS, fS,...
    delta0, omega0, e0, m0, ...
    v0, theta0, pg0, qg0, ...
    pref0, f0, ...
    Alpha,NetworkS)%LQRstep calculates the required feedback gain for the LQR step based on
%the desired steady-state zS and the previous steady-state z0. 
% The formulation solves an SDP per equation (18) but with zS known. 
% 
% Description of Outputs: 
% 1. K: the LQR gain, size(2*G, 3*G) 
% 2. trCostEstimate: is the estimate of the transient cost
% 
% Description of Inputs:
% 1. zS: the desired steady-state z=(x,a,u), size(2*N+7*G,1).
% 2. z0: the previous steady-state z0=(x0,a0,u0), size(2*N+7*G,1).
% 3. alpha

global OMEGAS Sbase N G L NodeSet GenSet LoadSet YMat GMat BMat Cg...
    YffVec YftVec  YtfVec YttVec

%  indices [these  do not change]
global deltaIdx omegaIdx eIdx mIdx  ...
    thetaIdx vIdx pgIdx qgIdx  fIdx prefIdx  SlackIdx

% machine [these do not change]
global  TauVec XdVec XqVec XprimeVec DVec MVec TchVec FreqRVec...

global Tlqr




%% Obtaining the jacobians:
 [ gx,ga,gu ] = gFunctionJacobVectorized(delta0, omega0, e0, m0, v0, theta0, pg0, qg0);
[ hx, ha] = hFunctionJacobVectorized( delta0, omega0, e0, m0, v0, theta0, pg0, qg0);
haInv=mldivide(ha,speye(size(ha)));
Asys=gx-ga*haInv*hx;
Bsys=gu;


%% costs:




x0=[delta0;omega0;e0;m0];



xs=[deltaS;omegaS;eS;mS];






cvx_begin quiet
cvx_solver sdpt3
variables   gama2
variable P(4*G,4*G) symmetric
variable Y(2*G,4*G) 
expression Qinv(4*G,4*G) 
expression Rinv(2*G,2*G) 
minimize( gama2) ;
subject to:







[Qinv,Rinv]=QinvRinv(pgS,qgS,Alpha, Qinv, Rinv,NetworkS);


-[-gama2, xs.'-x0.'; xs-x0, -P] ==semidefinite(4*G+1);

 -[Asys*P + P*Asys.' + Bsys*Y + Y.' * Bsys.', P , Y.';
    P, -Qinv, zeros(4*G, 2*G); 
    Y, zeros(2*G, 4*G), -Rinv] == semidefinite(10*G);

   
cvx_end

K=-Rinv*Bsys.'*conj(inv(P));
TrCostEstimate=(Tlqr/2)*gama2;

Gamma=(xs-x0).'*inv(P)*(xs-x0);



end


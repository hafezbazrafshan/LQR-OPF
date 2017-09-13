function [ gx,ga] = gTildeFunctionJacobVectorized(...
     delta, omega, e,m,...
     v,theta,pg,qg)
%GFUNCTIONJACOBVECTORIZED calculates the jacobian of g(x,a,u); 
% [gx,ga,gu]= gFunctionJacobVectorized(z, m_vec, d_vec, tau_vec, xd_vec,xprime_vec)
% calculates the jacobian of g(x,a,u). 
% 
% Description of Outputs: 
% 1. gx: the jacobian of g with respect to x, size(4*G,4*G); 
% 2. ga: the jacobian of g with respect to a, size(4*G, 2*N+2*G)
% 3. gu: the jacobian of g with respect to u, size(4*G, 2*G)
% 
% Description of Inputs: 
% 1. z: vector of z=(x,a,u) combining states, algebraic and control
% variables
% 2. m_vec: vector of generator inertia constants in pu, size(G,1). 
% 3. d_vec: vector of damping coefficients, size(G,1)
% 4. tau_vec: vector of direct axis transient open-circuit time constant,
% size(G,1). 
% 5. xd_vec:  vector of direct axis synchronous reactance pu, size(G,1).
% 6. xprime_vec: direct axis transient reactance pu, size(G,1).
% 7. Tch_vec
% 8. freqR_vec
% See also gFunctionJacob



global N G  GenSet LoadSet ...
    deltaIdx omegaIdx eIdx mIdx...
    thetaIdx vIdx pgIdx qgIdx prefIdx fIdx 

% LQR control
global KLQRstep

M=diag(MVec);
D=diag(DVec);
T=diag(TauVec);
Xprime=diag(XprimeVec); 
Xd=diag(XdVec);
Tch=diag(TchVec); 
freqR=diag(FreqRVec); 

InvM=diag(1./MVec);
InvT=diag(1./TauVec);
InvXprime=diag(1./XprimeVec);
InvTch=diag(1./TchVec);
InvfreqR=diag(1./FreqRVec);





Vg=v(GenSet);
V_g=diag(Vg);
thetag=theta(GenSet);
VgIdx=vIdx(GenSet);
thetagIdx=thetaIdx(GenSet);


[gx,ga]=









end


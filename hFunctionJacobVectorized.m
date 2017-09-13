function [ hx, ha ] = hFunctionJacobVectorized( delta, omega, e, m, v, theta, pg, qg )
%HFUNCTIONJACOBVECTORIZED calculates the jacobian of g(x,a,u); 
% [ hx, ha,hu ] = hFunctionJacob( z,...
 %    xq_vec, xprime_vec ) calculates the jacobian of h. 
% 
% Description of Outputs: 
% 1. hx: the jacobian of h with respect to x, size(2*N+2*G,4*G); 
% 2. ha: the jacobian of h with respect to a, size(2*N+2*G, 2*N+2*G)
% 3. hu: the jacobian of h with respect to u, size(2*N+2*G, 2*G)
% 
% Description of Inputs: 
% 1. z: vector of z=(x,a,u) combining states, algebraic and control
% variables, size(2*N+7*G,1).
% 5. xq_vec:  vector of quadrature axis synchronous reactance (pu) size(G,1)
% 6. xprime_vec: direct axis transient reactance pu, size(G,1).
%
% See also hFunctionJacob

global OMEGAS Sbase N G L NodeSet GenSet LoadSet YMat GMat BMat Cg...
    YffVec YftVec  YtfVec YttVec

global deltaIdx omegaIdx eIdx mIdx  ...
    thetaIdx vIdx pgIdx qgIdx  fIdx prefIdx  


% machine [these do not change]
global  TauVec XdVec XqVec XprimeVec DVec MVec TchVec FreqRVec...

% quantities:






Vg=v(GenSet);
thetag=theta(GenSet);


% constants
Xprime=diag(sparse(XprimeVec)); 
Xq=diag(sparse(XqVec));
XqInv=mldivide(Xq,speye(size(Xq)));
XprimeInv=mldivide(Xprime,speye(size(Xprime)));

% Matrix variables
E=diag(sparse(e));
V_g=diag(sparse(Vg));
VMat=diag(sparse(v));
CosMat=diag( cos(sparse(theta)));
SinMat=diag( sin(sparse(theta)));


% Indices
VgIdx=vIdx(GenSet);
thetagIdx=thetaIdx(GenSet);


h1Idx=(1:G).';
h2Idx=(h1Idx(end)+1:h1Idx(end)+G).';
h3Idx=(h2Idx(end)+1:h2Idx(end)+G).';
h4Idx=(h3Idx(end)+1:h3Idx(end)+G).';
h5Idx=(h4Idx(end)+1:h4Idx(end)+L).';
h6Idx=(h5Idx(end)+1:h5Idx(end)+L).';






hx=sparse(2*N+2*G, length([delta; omega;e; m])); 
ha=sparse(2*N+2*G,length([v;theta;pg;qg]));



hx(h1Idx,deltaIdx)= XprimeInv*E*V_g*diag(sparse(cos(delta-thetag)))+...
   XqInv*XprimeInv*(Xprime-Xq)*V_g*V_g*diag(sparse(cos(2*(delta-thetag))));
hx(h1Idx,eIdx)=XprimeInv*V_g*diag(sparse(sin(delta-thetag)));
hx(h2Idx,deltaIdx)=-XprimeInv*E*V_g*diag(sparse(sin(delta-thetag)))...
    -XqInv*XprimeInv*(Xprime-Xq)*V_g*V_g*diag(sparse(sin(2*(delta-thetag))));
hx(h2Idx,eIdx)=XprimeInv*V_g*diag(sparse(cos(delta-thetag)));



ha(h1Idx,pgIdx)=-speye(G); 
ha(h1Idx,VgIdx)= XprimeInv*E*diag(sparse(sin(delta-thetag)))+...
    XqInv*XprimeInv*(Xprime-Xq)*V_g*diag(sparse(sin(2*(delta-thetag))));
ha(h1Idx,thetagIdx)=-XprimeInv*E*V_g*diag(sparse(cos(delta-thetag)))-...
    XqInv*XprimeInv*(Xprime-Xq)*V_g*V_g*diag(sparse(cos(2*(delta-thetag))));

ha(h2Idx,qgIdx)=-speye(G); 
ha(h2Idx,VgIdx)=XprimeInv*E*diag(sparse(cos(delta-thetag)))- ...
    XqInv*XprimeInv*(Xprime+Xq)*V_g+...
   XqInv*XprimeInv*(Xprime-Xq)*V_g*diag(sparse(cos(2*(delta-thetag))));
ha(h2Idx,thetagIdx)=XprimeInv*E*V_g*diag(sparse(sin(delta-thetag)))+...
    XqInv*XprimeInv*(Xprime-Xq)*V_g*V_g*diag(sparse(sin(2*(delta-thetag))));


h35=sparse(N,length([v;theta;pg;qg]));
h35(:, pgIdx)= -Cg;
h35(:,vIdx)= diag(CosMat*GMat* CosMat*v) +VMat*CosMat*GMat*CosMat...
    - diag(CosMat*BMat*SinMat*v)- VMat*CosMat*BMat*SinMat...
    +diag( SinMat*BMat*CosMat*v)+ VMat*SinMat*BMat*CosMat...
    +diag(SinMat*GMat*SinMat*v)+ VMat*SinMat*GMat*SinMat;

h35(:,thetaIdx)=-diag(VMat*GMat*VMat*cos(theta))*SinMat- CosMat*VMat*GMat*VMat*SinMat...
    +diag(VMat*BMat*VMat*sin(theta))*SinMat-CosMat*VMat*BMat*VMat*CosMat...
    +diag(VMat*BMat*VMat*cos(theta))*CosMat-SinMat*VMat*BMat*VMat*SinMat...
    +diag(VMat*GMat*VMat*sin(theta))*CosMat+SinMat*VMat*GMat*VMat*CosMat;


ha(h3Idx,:)= h35(GenSet,:);
ha(h5Idx,:)=h35(LoadSet,:); 

h46=sparse(N,length([v;theta;pg;qg]));
h46(:, qgIdx)= -Cg;
h46(:,vIdx)= diag(SinMat*GMat*CosMat*v) +VMat*SinMat*GMat*diag(cos(theta))...
    - diag(SinMat*BMat*SinMat*v)- VMat*SinMat*BMat*diag(sin(theta))...
    - diag(CosMat*BMat*CosMat*v)- VMat*CosMat*BMat*diag(cos(theta))...
    -diag(CosMat*GMat*SinMat*v)-VMat*CosMat*GMat*diag(sin(theta));

h46(:,thetaIdx)=diag(VMat*GMat*VMat*cos(theta))*CosMat- SinMat*VMat*GMat*VMat*SinMat...
    -diag(VMat*BMat*VMat*sin(theta))*CosMat-SinMat*VMat*BMat*VMat*CosMat...
    +diag(VMat*BMat*VMat*cos(theta))*SinMat+CosMat*VMat*BMat*VMat*SinMat...
    +diag(VMat*GMat*VMat*sin(theta))*SinMat-CosMat*VMat*GMat*VMat*CosMat;
     
ha(h4Idx,:)=h46(GenSet,:);
ha(h6Idx,:)=h46(LoadSet,:);







end


function [ H,JH] =hAlgebraic(delta,omega, e,m, v,theta, pd, qd)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

% system constants [these do not change]
global OMEGAS Sbase N G L NodeSet GenSet LoadSet YMat GMat BMat Cg...
    YffVec YftVec  YtfVec YttVec

%  indices [these  do not change]
global deltaIdx omegaIdx eIdx mIdx  ...
    thetaIdx vIdx pgIdx qgIdx  fIdx prefIdx  SlackIdx

% machine [these do not change]
global  TauVec XdVec XqVec XprimeVec DVec MVec TchVec FreqRVec...
    

% dynamical simulations 
global TFinal TPert FSample NSamples NPertSamples Mass...
PertSet PPertValues QPertValues NoiseVarianceSet 



Xprime=diag(sparse(XprimeVec)); 
Xq=diag(sparse(XqVec));
XqInv=mldivide(Xq,speye(size(Xq)));
XprimeInv=mldivide(Xprime,speye(size(Xprime)));

vg=v(GenSet);
thetag=theta(GenSet);

EMat=diag(sparse(e));
V_g=diag(sparse(vg));
VMat=diag(sparse(v));
CosMat=diag( cos(sparse(theta)));
SinMat=diag( sin(sparse(theta)));



pg=XprimeInv*EMat*diag(sparse(sin(delta-thetag)))*vg+0.5*XqInv*XprimeInv*(Xprime-Xq)*diag(sparse(sin(2*(delta-thetag))))*V_g*vg;


% h35=VMat*CosMat*GMat*VMat*cos(theta)-VMat*CosMat*BMat*VMat*sin(theta)...
%   +VMat*SinMat*BMat*VMat*cos(theta)+VMat*SinMat*GMat*VMat*sin(theta)-Cg*pg;
h35=VMat*( CosMat*GMat*CosMat- CosMat*BMat*SinMat...
                    + SinMat*BMat*CosMat+SinMat*GMat*SinMat)*v- Cg*pg;
                
 qg=XprimeInv*EMat* diag(sparse(cos(delta-thetag)))*vg-...
     0.5*XqInv*XprimeInv*(Xprime+Xq)*V_g*vg+0.5*XqInv*XprimeInv*(Xprime-Xq)*diag(sparse(cos(2*(delta-thetag))))*V_g*vg;
 
 
h46=VMat*(SinMat*GMat*CosMat-SinMat*BMat*SinMat-CosMat*BMat*CosMat-CosMat*GMat*SinMat)*v-Cg*qg;                
%                h46=VMat*SinMat*GMat*VMat*cos(theta)-VMat*SinMat*BMat*VMat*sin(theta)...
%   -VMat*CosMat*BMat*VMat*cos(theta)-  VMat*CosMat*GMat*VMat*sin(theta)-Cg*qg;

h3=h35(GenSet)+pd(GenSet);
h5=h35(LoadSet)+pd(LoadSet);

h4=h46(GenSet)+qd(GenSet);
h6=h46(LoadSet)+qd(LoadSet);

H=[h3;h5;h4;h6];





Jpgvg=XprimeInv*EMat*diag(sparse(sin(delta-thetag)))+...
    XqInv*XprimeInv*(Xprime-Xq)*V_g*diag(sparse(sin(2*(delta-thetag))));
Jpgthetag=-XprimeInv*EMat*V_g*diag(sparse(cos(delta-thetag)))-...
    XqInv*XprimeInv*(Xprime-Xq)*V_g*V_g*diag(sparse(cos(2*(delta-thetag))));

Jqgvg=XprimeInv*EMat*diag(sparse(cos(delta-thetag)))- ...
    XqInv*XprimeInv*(Xprime+Xq)*V_g+...
   XqInv*XprimeInv*(Xprime-Xq)*V_g*diag(sparse(cos(2*(delta-thetag))));
Jqgthetag=XprimeInv*EMat*V_g*diag(sparse(sin(delta-thetag)))+...
    XqInv*XprimeInv*(Xprime-Xq)*V_g*V_g*diag(sparse(sin(2*(delta-thetag))));


Jh35=sparse(N,length([v;theta]));

Jh35(:,vIdx)= diag(CosMat*GMat* CosMat*v) +VMat*CosMat*GMat*CosMat...
    - diag(CosMat*BMat*SinMat*v)- VMat*CosMat*BMat*SinMat...
    +diag( SinMat*BMat*CosMat*v)+ VMat*SinMat*BMat*CosMat...
    +diag(SinMat*GMat*SinMat*v)+ VMat*SinMat*GMat*SinMat-Cg*Jpgvg*Cg.';
Jh35(:,thetaIdx)=-diag(VMat*GMat*VMat*cos(theta))*SinMat- CosMat*VMat*GMat*VMat*SinMat...
    +diag(VMat*BMat*VMat*sin(theta))*SinMat-CosMat*VMat*BMat*VMat*CosMat...
    +diag(VMat*BMat*VMat*cos(theta))*CosMat-SinMat*VMat*BMat*VMat*SinMat...
    +diag(VMat*GMat*VMat*sin(theta))*CosMat+SinMat*VMat*GMat*VMat*CosMat-Cg*Jpgthetag*Cg.';

Jh46=sparse(N,length([v;theta]));
 Jh46(:,vIdx)=diag(SinMat*GMat*CosMat*v) +VMat*SinMat*GMat*diag(cos(theta))...
    - diag(SinMat*BMat*SinMat*v)- VMat*SinMat*BMat*diag(sin(theta))...
    - diag(CosMat*BMat*CosMat*v)- VMat*CosMat*BMat*diag(cos(theta))...
    -diag(CosMat*GMat*SinMat*v)-VMat*CosMat*GMat*diag(sin(theta))-Cg*Jqgvg*Cg.';

Jh46(:,thetaIdx)=diag(VMat*GMat*VMat*cos(theta))*CosMat- SinMat*VMat*GMat*VMat*SinMat...
    -diag(VMat*BMat*VMat*sin(theta))*CosMat-SinMat*VMat*BMat*VMat*CosMat...
    +diag(VMat*BMat*VMat*cos(theta))*SinMat+CosMat*VMat*BMat*VMat*SinMat...
    +diag(VMat*GMat*VMat*sin(theta))*SinMat-CosMat*VMat*GMat*VMat*CosMat-Cg*Jqgthetag*Cg.';

Jh3=Jh35(GenSet,:); 
Jh5=Jh35(LoadSet,:); 
Jh4=Jh46(GenSet,:); 
Jh6=Jh46(LoadSet,:); 

JH=[Jh3;Jh5;Jh4;Jh6];


end


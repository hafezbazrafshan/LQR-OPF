function [ pFrom,pTo ] = determineLineFlows(lineIdx,v,theta)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here



% system constants [these do not change]
global OMEGAS Sbase N G L NodeSet GenSet LoadSet YMat GMat BMat Cg...
    YffVec YftVec  YtfVec YttVec

if isempty(lineIdx)
    pFrom=0;
   pTo=0;
   return;
end

fromBus=Network.branch(lineIdx,1); 
toBus=Network.branch(lineIdx,2); 
vFrom=v(fromBus); 
vTo=v(toBus); 

thetaFrom=theta(fromBus); 
thetaTo=theta(toBus); 


vFromComplex=vFrom.*exp(1j.*thetaFrom);
vToComplex=vTo.*exp(1j*thetaTo); 



iFrom=YffVec(lineIdx).*vFromComplex+YftVec(lineIdx).*vToComplex; 
iTo=YtfVec(lineIdx).*vFromComplex+YttVec(lineIdx).*vToComplex;


pFrom=real(vFromComplex.*conj(iFrom)); 
pTo=real(vToComplex.*conj(iTo)); 


end


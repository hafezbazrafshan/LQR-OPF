function [ pFrom,pTo ] = determineLineFlows(lineIdx,v,theta)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

global network yff_vec yft_vec ytf_vec ytt_vec


if isempty(lineIdx)
    pFrom=0;
   pTo=0;
   return;
end

fromBus=network.branch(lineIdx,1); 
toBus=network.branch(lineIdx,2); 
vFrom=v(fromBus); 
vTo=v(toBus); 

thetaFrom=theta(fromBus); 
thetaTo=theta(toBus); 


vFromComplex=vFrom.*exp(1j.*thetaFrom);
vToComplex=vTo.*exp(1j*thetaTo); 



iFrom=yff_vec(lineIdx).*vFromComplex+yft_vec(lineIdx).*vToComplex; 
iTo=ytf_vec(lineIdx).*vFromComplex+ytt_vec(lineIdx).*vToComplex;


pFrom=real(vFromComplex.*conj(iFrom)); 
pTo=real(vToComplex.*conj(iTo)); 


end


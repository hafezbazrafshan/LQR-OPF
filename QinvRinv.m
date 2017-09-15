function [ Qinv,Rinv ] = QinvRinv( pgS,qgS,alpha, Qinv, Rinv,NetworkS )
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here


global Sbase G omegaIdx deltaIdx eIdx mIdx OMEGAS








Qinv(sub2ind([4*G,4*G], omegaIdx, omegaIdx))=1*(1-alpha*pgS./(NetworkS.gen(:,9)./Sbase));
Qinv(sub2ind([4*G,4*G], deltaIdx, deltaIdx)) =1*(1-alpha*pgS./(NetworkS.gen(:,9)./Sbase));
Qinv(sub2ind([4*G, 4*G], eIdx,eIdx)) =1*(1-alpha*qgS./(NetworkS.gen(:,4)./Sbase));
Qinv(sub2ind([4*G,4*G],mIdx,mIdx))=1*(1-alpha*pgS./(NetworkS.gen(:,9)./Sbase));
Rinv(sub2ind([2*G, 2*G], 1:G,1:G))  =1*(1-alpha*pgS./(NetworkS.gen(:,9)./Sbase));
Rinv(sub2ind([2*G, 2*G], G+1:2*G, G+1:2*G)) =1*(1-alpha*qgS./(NetworkS.gen(:,4)./Sbase));

end


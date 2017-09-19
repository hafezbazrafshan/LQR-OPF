function [ Qinv,Rinv ] = QinvRinv( pgS,qgS,alpha, Qinv, Rinv,NetworkS )
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here


global Sbase G omegaIdx deltaIdx eIdx mIdx prefIdx fIdx OMEGAS GenSlackIdx





Qinv=eye(4*G); 
Rinv=eye(2*G);
% PIndices=setdiff(find(NetworkS.gen(:,9)>0),GenSlackIdx);
% QIndices=setdiff(find(NetworkS.gen(:,4)>0),GenSlackIdx);
PIndices=find(NetworkS.gen(:,9)>0);
QIndices=find(NetworkS.gen(:,4)>0);


Qinv(sub2ind([4*G,4*G], omegaIdx(PIndices), omegaIdx(PIndices)))=...
    1*(1.2-alpha*pgS(PIndices)./(NetworkS.gen(PIndices,9)./Sbase));
Qinv(sub2ind([4*G,4*G], deltaIdx(PIndices), deltaIdx(PIndices))) =...
    1*(1.2-alpha*pgS(PIndices)./(NetworkS.gen(PIndices,9)./Sbase));
Qinv(sub2ind([4*G, 4*G], eIdx(QIndices),eIdx(QIndices))) =...
    1*(1.2-alpha*qgS(QIndices)./(NetworkS.gen(QIndices,4)./Sbase));
Qinv(sub2ind([4*G,4*G],mIdx(PIndices),mIdx(PIndices)))=...
        1*(1.2-alpha*pgS(PIndices)./(NetworkS.gen(PIndices,9)./Sbase));
Rinv(sub2ind([2*G, 2*G], prefIdx(PIndices), prefIdx(PIndices)))  =...
        1*(1.2-alpha*pgS(PIndices)./(NetworkS.gen(PIndices,9)./Sbase));
Rinv(sub2ind([2*G, 2*G], fIdx(QIndices), fIdx(QIndices))) =...
        1*(1.2-alpha*qgS(QIndices)./(NetworkS.gen(QIndices,4)./Sbase));

% Qinv(sub2ind([4*G,4*G], omegaIdx(PIndices), omegaIdx(PIndices)))=...
%    1+ alpha*abs(pgS(PIndices)./(NetworkS.gen(PIndices,9)./Sbase));
% Qinv(sub2ind([4*G,4*G], deltaIdx(PIndices), deltaIdx(PIndices))) =...
%    1+alpha*abs(pgS(PIndices)./(NetworkS.gen(PIndices,9)./Sbase));
% Qinv(sub2ind([4*G, 4*G], eIdx(QIndices),eIdx(QIndices))) =...
%    1+alpha*abs(qgS(QIndices)./(NetworkS.gen(QIndices,4)./Sbase));
% Qinv(sub2ind([4*G,4*G],mIdx(PIndices),mIdx(PIndices)))=...
%       1+ alpha*abs(pgS(PIndices)./(NetworkS.gen(PIndices,9)./Sbase));
% Rinv(sub2ind([2*G, 2*G], prefIdx(PIndices), prefIdx(PIndices)))  =...
%        1+alpha*abs(pgS(PIndices)./(NetworkS.gen(PIndices,9)./Sbase));
% Rinv(sub2ind([2*G, 2*G], fIdx(QIndices), fIdx(QIndices))) =...
%        1+ alpha*abs(qgS(QIndices)./(NetworkS.gen(QIndices,4)./Sbase));


end


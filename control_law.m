function [pref,f]=control_law(ControlMode,delta,omega,e,m,...
                    v, theta, pg,qg,y)
      
                
           
% system constant variables
global G
    

    
    


    






% new equilibrium 
global deltaS omegaS eS mS  prefS fS...


% LQR control
global KLQRstep





                
     switch ControlMode
         
         
         case 'OpenLoop'
             
             
             pref=prefS;
             f=fS;
             
         case 'LQR'
             
             u=KLQRstep*[delta- deltaS; omega- omegaS; e-eS;m-mS]+[prefS;fS];
             
             
             
             pref=u(1:G); 
            f=u(G+1:end);


    


             case 'AGC'
              u=KLQRstep*[delta- deltaS; omega- omegaS; e-eS;m-mS]+[prefS;fS];
%               f=KLQRstep(G+1:end,eIdx)*[e-eS]+fS;
f=u(G+1:end);
%              pref=zeros(G,1);
             for ii=1:NumberOfAreas
%                  pref(GensPerArea{ii,1})=pref0(GensPerArea{ii,1})+0.5*(ParticipationFactors{ii}.*y(ii))+0.5*(pg(GensPerArea{ii,1})-pg0(GensPerArea{ii,1}));
%   pref(GensPerArea{ii,1})=pref0(GensPerArea{ii,1})+1*(ParticipationFactors{ii}.*y(ii))+0*(pg(GensPerArea{ii,1})-pg0(GensPerArea{ii,1}));
%   pref(GensPerArea{ii,1})=(ParticipationFactors{ii}.*y(ii))+0.1*(pg(GensPerArea{ii,1})-pgS(GensPerArea{ii,1}));
%   pref(GensPerArea{ii,1})=ParticipationFactors{ii}.*y(ii);
pref(GensPerArea{ii,1})=ParticipationFactors{ii}.*y(GensPerArea{ii,1});





             end

end


             
             


             
             
     
     

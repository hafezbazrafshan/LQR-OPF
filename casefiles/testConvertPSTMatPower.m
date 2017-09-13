

d2asbegp;
clearvars -except bus line mac_con
clc;

% determine the number of buses
N=size(bus,1); 

MatBus=zeros(N,13); 

% % % 
% % %  Bus Data Format
% % %        1   bus number (positive integer)
% % %        2   bus type
% % %                PQ bus          = 1
% % %                PV bus          = 2
% % %                reference bus   = 3
% % %                isolated bus    = 4
% % %        3   Pd, real power demand (MW)
% % %        4   Qd, reactive power demand (MVAr)
% % %        5   Gs, shunt conductance (MW demanded at V = 1.0 p.u.)
% % %        6   Bs, shunt susceptance (MVAr injected at V = 1.0 p.u.)
% % %        7   area number, (positive integer)
% % %        8   Vm, voltage magnitude (p.u.)
% % %        9   Va, voltage angle (degrees)
% % %        10  baseKV, base voltage (kV)
% % %        11  zone, loss zone (positive integer)
% % %        12  maxVm, maximum voltage magnitude (p.u.)
% % %        13  minVm, minimum voltage magnitude (p.u.)

Sbase=100;
for ii=1:N
    % determine column 1 (bus number)
    MatBus(ii,1)=round(bus(ii,1)); % bus number are the same
    
    % determine column 2 (bus type)
   if bus(ii,10)==1
       MatBus(ii,2)=3; % 3 indicates slack for MATPOWER
   elseif bus(ii,10)==2
       MatBus(ii,2)=2; % 2 indicates PV bus for MATPOWER
   elseif bus(ii,10)==3
       MatBus(ii,2)=1; % 1 indicates PQ bus for MATPOWER
   end
   
   
   % determine column 3 (real power demand)
   MatBus(ii,3)=bus(ii,6)*Sbase; 
   
   % determine column 4 (reactive power demand) 
   MatBus(ii,4)=bus(ii,7)*Sbase;
   
   % determine column 5 (shunt conductance) 
   MatBus(ii,5)=bus(ii,8)*Sbase;
   
   % determine column 6 (shunt susceptance) 
   MatBus(ii,6)=bus(ii,9)*Sbase;
   
   % determine column 7 (area) (not specified)
   MatBus(ii,7)=1;
   
   % determine column 8 (voltage magnitude)
   MatBus(ii,8)=bus(ii,2);
   
   % determine column 9 (voltage angle in degrees)
   MatBus(ii,9)=bus(ii,3);
   
   % determine column 10 (base KV)
   MatBus(ii,10)=500; % based on MATPOWER
   
   % determine column 11 (zone)
   MatBus(ii,11)=1;
   
   % determine column 12 (Vmax)
   MatBus(ii,12)=max(bus(:,2))+0.2;
   
   % determine column 13 (Vmin) 
   MatBus(ii,13)=min(bus(:,2))-0.2;
   
    
end



%    Generator Data Format
%        1   bus number
%    (-)     (machine identifier, 0-9, A-Z)
%        2   Pg, real power output (MW)
%        3   Qg, reactive power output (MVAr)
%        4   Qmax, maximum reactive power output (MVAr)
%        5   Qmin, minimum reactive power output (MVAr)
%        6   Vg, voltage magnitude setpoint (p.u.)
%    (-)     (remote controlled bus index)
%        7   mBase, total MVA base of this machine, defaults to baseMVA
%    (-)     (machine impedance, p.u. on mBase)
%    (-)     (step up transformer impedance, p.u. on mBase)
%    (-)     (step up transformer off nominal turns ratio)
%        8   status,  >  0 - machine in service
%                     <= 0 - machine out of service
%    (-)     (% of total VAr's to come from this gen in order to hold V at
%                remote bus controlled by several generators)
%        9   Pmax, maximum real power output (MW)
%        10  Pmin, minimum real power output (MW)
%    (2) 11  Pc1, lower real power output of PQ capability curve (MW)
%    (2) 12  Pc2, upper real power output of PQ capability curve (MW)
%    (2) 13  Qc1min, minimum reactive power output at Pc1 (MVAr)
%    (2) 14  Qc1max, maximum reactive power output at Pc1 (MVAr)
%    (2) 15  Qc2min, minimum reactive power output at Pc2 (MVAr)
%    (2) 16  Qc2max, maximum reactive power output at Pc2 (MVAr)
%    (2) 17  ramp rate for load following/AGC (MW/min)
%    (2) 18  ramp rate for 10 minute reserves (MW)
%    (2) 19  ramp rate for 30 minute reserves (MW)
%    (2) 20  ramp rate for reactive power (2 sec timescale) (MVAr/min)
%    (2) 21  APF, area participation factor

G=length(mac_con(:,1));
MatGen=zeros(G,21);

for ii=1:G
    % determine column 1:
    MatGen(ii,1)=mac_con(ii,2);
    
    % determine column 2 (real power generation)
    BusNumber=MatGen(ii,1); % the bus number
    MatGen(ii,2)=bus(BusNumber,4)*Sbase;
   
    % determine column 3 (reactive power generation)
    MatGen(ii,3)=bus(BusNumber,5).*Sbase;

    % determine column 6 (voltage magnitude)
    MatGen(ii,6)=bus(BusNumber,2); 
    
    
    
end

    
    % determine column 4,5 (max,min reactive power)
    MatGen(:,4)= max(bus(:,5))*Sbase+50;
    MatGen(:,5)=min(bus(:,5))*Sbase-50;
    
    % determine column 9,10 (max,min active power)
    MatGen(:,9)=max(bus(:,4))*Sbase+50; 
    MatGen(:,10)=0;
    
    % determine column 7
    MatGen(:,7)=Sbase;
    
    % determine column 8
    MatGen(:,8)=1;
    
    
%     Branch Data Format
%        1   f, from bus number
%        2   t, to bus number
%    (-)     (circuit identifier)
%        3   r, resistance (p.u.)
%        4   x, reactance (p.u.)
%        5   b, total line charging susceptance (p.u.)
%        6   rateA, MVA rating A (long term rating)
%        7   rateB, MVA rating B (short term rating)
%        8   rateC, MVA rating C (emergency rating)
%        9   ratio, transformer off nominal turns ratio ( = 0 for lines )
%            (taps at 'from' bus, impedance at 'to' bus,
%             i.e. if r = x = 0, then ratio = Vf / Vt)
%        10  angle, transformer phase shift angle (degrees), positive => delay
%    (-)     (Gf, shunt conductance at from bus p.u.)
%    (-)     (Bf, shunt susceptance at from bus p.u.)
%    (-)     (Gt, shunt conductance at to bus p.u.)
%    (-)     (Bt, shunt susceptance at to bus p.u.)
%        11  initial branch status, 1 - in service, 0 - out of service
%    (2) 12  minimum angle difference, angle(Vf) - angle(Vt) (degrees)
%    (2) 13  maximum angle difference, angle(Vf) - angle(Vt) (degrees)
%            (The voltage angle difference is taken to be unbounded below
%             if ANGMIN < -360 and unbounded above if ANGMAX > 360.
%             If both parameters are zero, it is unconstrained.)

L=size(line,1);
MatBranch=zeros(L,13);

MatBranch(:,13)=360;
MatBranch(:,12)=-360;
MatBranch(:,11)=1;

MatBranch(:,1)=round(line(:,1)); % from bus
MatBranch(:,2)=round(line(:,2)); % to bus
MatBranch(:,3)=line(:,3); % r
MatBranch(:,4)=line(:,4); % x
MatBranch(:,5)=line(:,5); % b

MatBranch(:,9)=line(:,6); % tap ratio
MatBranch(:,10)=line(:,7); % tap angle


mpc.bus=MatBus;
mpc.branch=MatBranch;
mpc.gen=MatGen;
mpc.baseMVA=Sbase;
mpc.gencost = repmat([
	2	0	0	3	0.01	0.3	0.2],G,1); 

savecase('CaseD2AsbegpFromPST',mpc);
    
    

[fb,ISI,Var,minavg,fbLT,ISILT,burnumLT,swnumLT]=TCcelloutput;

function [fb,ISI,Var,minavg,fbLT,ISILT,burnumLT,swnumLT] = TCcelloutput


t_int=[0,3000];

%%Low-ThreshBurst Data
burnum=16;
dc=.3;
fb=.0085;
swnum=2;
ISI=12;
Var=0;
minavg=-55;

LTBD=[fb,dc,burnum,swnum,ISI,Var,minavg];


% % % % Call parameterization function
[parameters,locmins] = neuron_parameterization(LTBD,t_int);

disp('Parameters')
disp(parameters)



[fb,dc,burnum,swnum,ISI,Var,minavg,timeSIM,VoltageSIM]=landscape(t_int,parameters);
[fbLT,dcLT,burnumLT,swnumLT,ISILT,VarLT,minavgLT,timeSIMLT,VoltageSimLT]=LowThresh(t_int,parameters);


figure(1)
plot(timeSIM,VoltageSIM)

figure(2)
plot(timeSIMLT,VoltageSimLT)

end



function [parameters,locmins] = neuron_parameterization(LTBD,t_int)


%[Iapp,gTLT,gcal,gM,gH,gKL,gK,gNAP,gL,gNA]

 % % Establish initial conditions for the ODE
%parms=[.1,4,.13,.04,.001,.001,4,.04,.01];% % Initial guess for parameters
%parms=[.27,4,.15,.04,.03,.001,4,.03,.01];


lb = [-.5,0,.1,0,0,0,0,.2,0,0];   % % lower bound vector
ub = [2,  4,.6,2,1,1,1,.8,1,1];   % % upper bound vector

A=zeros(1,10);
A(1,1)=1;
A(1,2)=-1;

B=0;
%parms=rand(1,length(ub)).*ub
%parms=[.03,.11,4.1,.14,.041,.0011,.0011,4.2,.041,.011];

% % % % fmincon is the optimization function used to produce MLE estimates
% % % % in our MLE procedure. The error function called is based on the NLL

options=gaoptimset("Display","iter",'OutputFcn',@gaoutfun,'CrossoverFraction',0.2);
[xga,fval,efga,outga,population,scores] = ga(@(x)error(LTBD,t_int,...
    [x(1) x(2) x(3) x(4) x(5) x(6) x(7) x(8) x(9) x(10)]),10,A,B,[],[],lb,ub,[],options);


 locmins=horzcat(population,scores);
toc

[M,I] = min(scores);




parameters=population(I,:);

end


function err=error(LTBD,t_int,parms)

[fb,dc,burnum,swnum,ISI,Var,minavg,time,V]=landscape(t_int,parms);
[fbLT,dcLT,burnumLT,swnumLT,ISILT,VarLT,minavgLT,timeSIMLT,VoltageSimLT]=LowThresh(t_int,parms);

ErrW=[100,5,20,1,20,1,2];
err1=ErrW(1)*((10-fb*1000)^2*(fb*1000>=8 && fb*1000<=12)+(20*(10-fb*1000)^2-76)*(fb*1000<8 || fb*1000>12));
err2=ErrW(2)*((16-ISI)^2*(ISI>=14 && ISI<=18)+(20*(16-ISI)^2-76)*(ISI<14 || ISI>18));
err3=ErrW(3)*((55+minavg)^2*(minavg>=-57 && minavg<=-53)+(20*(55+minavg)^2-76)*(minavg<-57 || minavg>-53));
err4=ErrW(4)*Var;

err5=ErrW(5)*((2-fbLT*1000)^2*(fbLT*1000>=1 && fbLT*1000<=3)+(20*(2-fbLT*1000)^2-19)*(fbLT*1000<1 || fbLT*1000>3));
err6=ErrW(6)*((7-ISILT)^2*(ISILT>=5 && ISILT<=9)+(20*(7-ISILT)^2-76)*(ISILT<5 || ISILT>9));
err7=ErrW(7)*(burnumLT-swnumLT)^2;
err=err1+err2+err3+err4+err5+err6;


% ErrW=[100,5,10,1,20,1,2];
% err1=ErrW(1)*((7-fb*1000)^2*(fb*1000>=5 && fb*1000<=9)+(20*(7-fb*1000)^2-76)*(fb*1000<5 || fb*1000>9));
% err2=ErrW(2)*((16-ISI)^2*(ISI>=14 && ISI<=18)+(20*(16-ISI)^2-76)*(ISI<14 || ISI>18));
% err3=ErrW(3)*((50+minavg)^2*(minavg>=-55 && minavg<=-45)+(20*(50+minavg)^2-475)*(minavg<-55 || minavg>-45));
% err4=ErrW(4)*Var;
% 
% err5=ErrW(5)*(2*(1.5-fbLT*1000)^2*(fbLT*1000<=1.5)+(-20*(1.5-fbLT*1000)^3)*(fbLT*1000>1.5));
% err6=ErrW(6)*((7-ISILT)^2*(ISILT>=5 && ISILT<=9)+(20*(7-ISILT)^2-76)*(ISILT<5 || ISILT>9));
% err7=ErrW(7)*(burnumLT-swnumLT)^2;
% err=err1+err2+err3+err4+err5+err6;
%if i do run to cursors on err, it spits out the parameters that it has
%.1*(Var-LTBD(6))^2

end

function [fb,dc,burnum,swnum,ISI,Var,minavg2,timeSIM,VoltageSim]=landscape(t_int,parms)

Cond=2;
[time,S]=voltageclamp(t_int,parms,Cond);

timeSIM=time;
VoltageSim=S(:,1)';

Tsp=-20;
%V=[S(:,1)',NaN]';

TIMESS=time;
runtoSS = time < 1000;
TIMESS(runtoSS)=[];
time=TIMESS;
% 
V=S(:,1)';
V(runtoSS)=[];
% 
% Tsp=-20;
V=[V,NaN]';
Spikes=[];
sw=[];
tsw=-65;


for i= 1:length(time)
   if V(i)<=Tsp && V(i+1)>Tsp
       Spikes=[Spikes,time(i)];
   end
   if V(i)>=tsw && V(i+1)<tsw
       sw=[sw,time(i)];
   end   
end
swnum=length(sw);


deltaspt=35;
Spikes=[time(1),Spikes,time(end)];
sdp=diff(Spikes(2:end));
sdm=diff(Spikes(1:end-1));

x=sdp;
bigNumbersLocations = x > deltaspt;
x(bigNumbersLocations) = [];
ISI=mean(x);
Var=var(x);
%do something with the variance, make the variance tightly knit

bs=[];
be=[];

for k=1:length(sdp)
    if sdp(k)<deltaspt && sdm(k) >deltaspt
        bs=[bs,Spikes(k+1)];
    elseif sdp(k) <deltaspt && Spikes(k)==time(1)
        bs=[bs,Spikes(k+1)];         
    end
    if sdp(k)>deltaspt && sdm(k) <deltaspt
       be=[be,Spikes(k+1)];
       
    elseif sdm(k) <deltaspt && Spikes(k+1)==Spikes(end-1)
       be=[be,Spikes(k+1)]; 
    end
end

%disp(Spikes)
% disp(sdm)
% %disp(sdp)
%   disp(bs);
%   disp(be);

if length(be)==length(bs)-1
     bs(end)=[];
end

if length(bs)==length(be)-1
     be(1)=[];
end
%  disp(bs);
%  disp(be);

if length(bs)~=length(be)
    disp(parms)
    disp(bs)
    disp(be)
end


delb=be-bs;
taub=diff(bs);
if length(bs)==1
    fb=0;
    dc=0;
else
    fb=mean(1./taub);
    dc=mean(delb)/mean(taub);
end

burnum=length(be);

minIDs = islocalmin(V);
dummy1=time(minIDs);
dummy2=V(minIDs);
TF=islocalmin(dummy2);
%display(dummy2(TF));
minavg=mean(dummy2(TF));

M=sort(dummy2(TF));

if burnum>1 && length(M)>=burnum
    minavg2=mean(M(1:(burnum-1)));
else
    minavg2=mean(dummy2(TF));
end

end

function [fbLT,dcLT,burnumLT,swnumLT,ISILT,VarLT,minavgLT,timeSIMLT,VoltageSimLT]=LowThresh(t_int,parms)

Cond=1;
%[time,S]=voltageclamp();
[time,S]=voltageclamp(t_int,parms,Cond);

timeSIMLT=time;
VoltageSimLT=S(:,1)';

Tsp=-20;
%V=[S(:,1)',NaN]';

TIMESS=time;
runtoSS = time < 1000;
TIMESS(runtoSS)=[];
time=TIMESS;
% 
V=S(:,1)';
V(runtoSS)=[];
% 
% Tsp=-20;
V=[V,NaN]';
Spikes=[];
sw=[];
tsw=-65;


for i= 1:length(time)
   if V(i)<=Tsp && V(i+1)>Tsp
       Spikes=[Spikes,time(i)];
   end
   if V(i)>=tsw && V(i+1)<tsw
       sw=[sw,time(i)];
   end   
end
swnumLT=length(sw);


deltaspt=50;
Spikes=[time(1),Spikes,time(end)];
sdp=diff(Spikes(2:end));
sdm=diff(Spikes(1:end-1));

x=sdp;
bigNumbersLocations = x > deltaspt;
x(bigNumbersLocations) = [];
ISILT=mean(x);
VarLT=var(x);
%do something with the variance, make the variance tightly knit

bs=[];
be=[];

for k=1:length(sdp)
    if sdp(k)<deltaspt && sdm(k) >deltaspt
        bs=[bs,Spikes(k+1)];
    elseif sdp(k) <deltaspt && Spikes(k)==time(1)
        bs=[bs,Spikes(k+1)];         
    end
    if sdp(k)>deltaspt && sdm(k) <deltaspt
       be=[be,Spikes(k+1)];
       
    elseif sdm(k) <deltaspt && Spikes(k+1)==Spikes(end-1)
       be=[be,Spikes(k+1)]; 
    end
end

%disp(Spikes)
% disp(sdm)
% %disp(sdp)
%   disp(bs);
%   disp(be);

if length(be)==length(bs)-1
     bs(end)=[];
end

if length(bs)==length(be)-1
     be(1)=[];
end
%  disp(bs);
%  disp(be);
if length(bs)~=length(be)
    disp(parms)
    disp(bs)
    disp(be)
end


delb=be-bs;
taub=diff(bs);
if length(bs)==1
    fbLT=0;
    dcLT=0;
else
    fbLT=mean(1./taub);
    dcLT=mean(delb)/mean(taub);
end

burnumLT=length(be);

minIDs = islocalmin(V);
dummy2=V(minIDs);
TF=islocalmin(dummy2);
%display(dummy2(TF));
minavgLT=mean(dummy2(TF));

end

function [time,S]=voltageclamp(t,parms,Cond)

S0 =[-80,.05,.6,.3,.0004,.3,.3,.05,.6,.05,.6]; 
options = odeset('RelTol',1e-6,'AbsTol',1e-6); 
[time,S]=ode15s(@(t,S)deRHS(t,S,parms,Cond),t,S0,options);
S(end,:);

end

function s_prime=deRHS(t,S,parms,Cond) 

%%%FITTED PARMS
IappLT=parms(1);
IappHT=parms(2);
gTLT_IN=parms(3);
gcal=parms(4);
gM=parms(5);
gH_IN =parms(6);
gKL_IN=parms(7);
gK_IN=parms(8);
gNAP=parms(9);
gL_IN=parms(10);
%gNa_IN=parms(11);

    if Cond==1
        Iapp=IappLT;
    elseif Cond==2
        Iapp=IappHT;
    end

%Fixed Parms
gNa_IN=90;
TEMP = 310.35;
vshift=55;
F=96485.3;
Rgas=8.314;
z = 2;
EH_IN = -43;
EKL_IN=-100;
EL_IN=-85;
EK_IN=-100;
ENa_IN=50;
Cm_IN = 1 ;

%State Vars
V=S(1);
mNa=S(2);
hNa=S(3);
mK=S(4);
Cai=S(5);
mTLT=S(6);
hTLT=S(7);
mcal=S(8);
hcal=S(9);
r=S(10);
p=S(11);

%Rates
mV=1;
ms=1;
Vt=V+vshift;
alphamNa=(0.32*(13*mV-Vt))./(exp((13*mV-Vt)/4/mV)-1);
betamNa = (0.28*(Vt-40*mV))./(exp((Vt-40*mV)/5/mV)-1);
alphahNa=(0.128*exp((17*mV-Vt)/18)) ;
betahNa = 4./(1+exp((40*mV-Vt)/5));
mNainf = alphamNa./(alphamNa+betamNa) ;
taumNa=1*ms./(alphamNa+betamNa)  ;   
hNainf = alphahNa./(alphahNa+betahNa); 
tauhNa=1*ms./(alphahNa+betahNa); 
alphamK = (0.032*(15*mV-Vt)/mV)./(exp((15*mV-Vt)/5/mV)-1);
betamK=0.5 * exp((10*mV-Vt)/40/mV);
mKinf = alphamK./(alphamK+betamK);
taumK=1*ms/(alphamK+betamK);
tausr=(1*ms)./(exp(-.086*V-14.59)+exp((.07*V-1.87)));
rinf = 1./(1+exp((V+75*mV)/5.5));
ECA_IN=Rgas*TEMP/(z*F)*log(2/Cai)*1000;


hTLTinf = 1/(1+exp((V+86*mV)/4/mV));
%tauhTLT=  1*ms*exp((V+470*mV)/66.6/mV).*heaviside(-V-83)+1*ms*(exp(-(V+25*mV)/10.5)+28).*heaviside(V+83);
taumTLT= .612*ms+(1*ms)./(exp(-(V+135)/16.7)+(exp((V+19.8)/18.2)));
mTLTinf = 1./(1+exp(-(V+62*mV)/6.2/mV));

pinf=1./(1+exp(-(V+35)/10));
taup=1000./(3.3*exp((V+35)/20)+exp(-(V+35)/20));

mNAPinf=1./(1+exp(-(V+50)/5));

if V <-83
    tauhTLT=1*ms*exp((V+470*mV)/66.6/mV);
else
    tauhTLT=1*ms*(exp(-(V+25*mV)/10.5)+28);
end

alphaq=.055*(-27-V)./(exp((-27-V)/3.8)-1);
betaq=.94*exp((-75-V)/17);
alphar=.000457*exp((-13-V)/50);
betar=.0065/(exp((-15-V)/28)+1);


%Currents
INa=gNa_IN*mNa.^3.*hNa.*(V-ENa_IN);
IK=10*gK_IN.*mK.^4.*(V-EK_IN);
IL=.02*gL_IN*(V-EL_IN);
IKL=.02*gKL_IN*(V-EKL_IN);
IH = .3*gH_IN * r .* (V-EH_IN);
ITLT = 10*gTLT_IN * mTLT.^2 .* hTLT .* (V-ECA_IN);
ICAL= gcal*mcal.^2.*hcal.*(V-ECA_IN);
IM=.3*gM*p.*(V-EK_IN);
INAP=.5*gNAP*mNAPinf.*(V-ENa_IN);

%%%RHS of the ODE equations

%-INa-IK-IL-IKL-IH-ITLT+Iapp-ITHT-ICAL-IAHP-ICAN

s_prime=[(-INa-IK-IL-IKL-IH-ITLT+Iapp-ICAL-IM-INAP)/Cm_IN,...
    (mNainf-mNa)./taumNa,...
    (hNainf-hNa)./tauhNa,...
    (mKinf-mK)./taumK,...
    ((-ITLT-ICAL)/(2*96485.3*.5))- (Cai-0.00005)/10,...
    4.6*(mTLTinf-mTLT)./taumTLT,...
    3.7*(hTLTinf-hTLT)./tauhTLT,...
    alphaq.*(1-mcal)-betaq.*mcal,...
    alphar.*(1-hcal)-betar.*hcal,...
    (rinf-r)./tausr,...
    (pinf-p)./taup]';
  
  
end

function [state,options,optchanged] = gaoutfun(options,state,flag)
persistent history 
optchanged = false;
    switch flag
    case 'init'
        M=state.Score;
        [sortedM,linearIndex] = sort(M(:));
        [rowIndex,~] = ind2sub(size(M),linearIndex);
        n = 15;
        sortedM(1:n);
        toprows=rowIndex(1:n);
%         Rows=state.Population(toprows,:);
%         disp([state.Population state.Score])
%         disp([state.Population(toprows,:) state.Score(toprows,:)])
        %BESTROWS=[state.Population(toprows,:) state.Score(toprows,:)];
        BESTROWS=[state.Population(toprows,:)];
        
        history(:,:,1) = BESTROWS;
        assignin('base','gapopulationhistory',history);
    case 'iter'
        % Update the history every 10 generations.
        if rem(state.Generation,3) == 0
            M=state.Score;
            [sortedM,linearIndex] = sort(M(:));
            [rowIndex,~] = ind2sub(size(M),linearIndex);
            n = 15;
            sortedM(1:n);
            toprows=rowIndex(1:n);
            BESTROWS=[state.Population(toprows,:)];
            
            ss = size(history,3);
            history(:,:,ss+1) = BESTROWS;
            assignin('base','gapopulationhistory',history);
        end
        % Find the best objective function, and stop if it is low.
        
        
        % Update the fraction of mutation and crossover after 25 generations.
         if state.Generation == 40
             options.CrossoverFraction = 0.8;
             optchanged = true;
         end
    case 'done'
        % Include the final population in the history.
        M=state.Score;
            [sortedM,linearIndex] = sort(M(:));
            [rowIndex,~] = ind2sub(size(M),linearIndex);
            n = 15;
            sortedM(1:n);
            toprows=rowIndex(1:n);
            BESTROWS=[state.Population(toprows,:)];
            
            ss = size(history,3);
            history(:,:,ss+1) = BESTROWS;
            assignin('base','gapopulationhistory',history);
    end
end
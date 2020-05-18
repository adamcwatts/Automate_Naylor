% ConvNaylor
% convolution trial
% import time related data
% Naylor Data Exel file
Time=Time2;
HeatFlux=Heatflux1;
RH=RH2;
Temp=Temp2;
% change these 4 to the specific columns and worksheets you want to analyze
% manual - not automated


HF0= mean(HeatFlux(1:20));
ff= 1.0;     % if ff=1, only 1 exotherm timeconstant used in model
tau1= 25;
tau2= 5000;  %  2nd time constant if needed
K= 4.78;
tauReg= 35;   
fs= 1-ff;    % fraction of 2nd time constant if needed

tti = length(Time);
Tend = Time(tti-1);
deltaTime= Tend/length(Time);

RHdiff(1)=0;

Naylor = 1*(ff*(1-exp(-Time/tau1)) + fs*(1-exp(-Time/tau2))).*K.*exp(-Time/tauReg);
for i = 3:tti,
    RHdiff(i) = RH(i)-RH(i-1);
end
tester = conv(RHdiff,Naylor,'full');
nn=length(tester);
convTime = 0:Tend/nn:Tend;
ResVec = interp1(convTime(1:nn/2),tester(1:nn/2),Time);
initTemp = mean(Temp(1:20))
TempEffect = HF0*initTemp./Temp;
HFmodel = TempEffect - ResVec;
figure(1);
plot(Time, HeatFlux,'r',Time,RH,'b',Time,10*Naylor,'r',2*Time,HFmodel,'k',Time,6*Temp,'g');

deltaRH = RH(tti-1)-RH(1);
ExothermicEffect = deltaRH*sum(Naylor(1:tti-1))*deltaTime
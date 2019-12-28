%% JT-60SA closed_loop_simulator
%clear all
close all


distFlag = 1;
stabFlag = 1;
startTime = 0;
stopTime = 10; 
switchTime = startTime-0.01;
cdTime = stopTime; % no current drive
neg_time=2;
transitionTime = 1.5;
interpTime = [startTime - neg_time; startTime; startTime+transitionTime; stopTime];
%interpTimeW = [startTime - neg_time; startTime; startTime+transitionTime-1; startTime+transitionTime; stopTime];

nPF = 10;

% Number of gaps
nOfGaps = 85;
% Number of strike point
nOfStrikes = 8;
% Number of flux sensors
nOfFluxSens = 34;
% Number of magnetic sensors
nOfMagneticSens = 45;
%%% Select number of control points
% Number of controlled fluxes
nOfFluxCntrl = nOfFluxCntrl +1;
% nOfFluxCntrl = 20;
% nOfFluxCntrl = 9;
%  nOfFluxCntrl = 7;
 % :D

% Load the initial equilibrium (which is also the linear model that will
% be used for the simulation

     
 if nOfFluxCntrl == 9
    simConf = '../../L-models_90deg/SOF@18d66s_japanese_sensors_8cntrlFluxPnts_CL.mat';
    load('../../L-models_90deg/SOF@18d66s_japanese_sensors_8cntrlFluxPnts.mat','Input_struct','x_np','y_np','y_type')
    targetConf = '../../L-models_90deg/SOF@18d66s_japanese_8cntrlFluxPntsSqzd_CL.mat';
 targetEquil = load('../../L-models_90deg/SOF@18d66s_japanese_8cntrlFluxPntsSqzd.mat');
 
 elseif nOfFluxCntrl ==7
      simConf = '../../L-models_90deg/SOF@18d66s_japanese_sensors_6cntrlFluxPnts_CL.mat';
      load('../../L-models_90deg/SOF@18d66s_japanese_sensors_6cntrlFluxPnts.mat','Input_struct','x_np','y_np','y_type')
     targetConf = '../../L-models_90deg/SOF@18d66s_japanese_6cntrlFluxPntsSqzd_CL.mat';
     targetEquil = load('../../L-models_90deg/SOF@18d66s_japanese_6cntrlFluxPntsSqzd.mat');

     
 else %%% 19 gaps
     
         simConf = '../../L-models_90deg/SOF@18d66s_japanese_sensors_cntrlFluxPnts_CL.mat';
        load('../../L-models_90deg/SOF@18d66s_japanese_sensors_cntrlFluxPnts.mat','Input_struct','x_np','y_np','y_type')
        targetConf = '../../L-models_90deg/SOF@18d66s_japanese_cntrlFluxPntsSqzd_2_CL.mat';
        targetEquil = load('../../L-models_90deg/SOF@18d66s_japanese_cntrlFluxPntsSqzd_2.mat');
 end

load(simConf);
targetConf = load(targetConf);



LinearModel.R(LinearModel.PlasmaCurrentInfo.StatePosition, LinearModel.PlasmaCurrentInfo.StatePosition) = 0; % force plasma resistance to 0

%% Initialize gaps
gapNames = {};
strikeNames = {};
FluxSensNames= {};
MagnSensNames={};

% Gap and strike names
for i = 1:nOfGaps
    gapNames{i} = sprintf('GAP%02d',i);
end

for i = 1:nOfStrikes
    strikeNames{i} = sprintf('GAP%02d',nOfGaps+i);
end

for i = 1:nOfFluxSens
    FluxSensNames{i} = sprintf('Flux_%03d',i);
end

for i = 1:nOfMagneticSens
    MagnSensNames{i} = sprintf('Bpol_%03d',i);
end


if nOfFluxCntrl == 9
    for i = 1:nOfFluxCntrl
     FluxCntrlNames{i} = sprintf('FluxCntrl_%03d',i);
     end
 elseif nOfFluxCntrl ==7
      for i = 1:nOfFluxCntrl
    FluxCntrlNames{i} = sprintf('FluxCntrl_%03d',i);
     end
else
     for i = 1:nOfFluxCntrl
    FluxCntrlNames{i} = sprintf('FluxCntrl_%03d',i+34);
end
end





% Gap and strike and sensors indexes
gapIdx = signalIndexByName(gapNames,LinearModel.OutputsInfo.Name,LinearModel.OutputsInfo.OutputPosition);
strikeIdx = signalIndexByName(strikeNames,LinearModel.OutputsInfo.Name,LinearModel.OutputsInfo.OutputPosition);
FluxSensIdx = signalIndexByName(FluxSensNames,LinearModel.OutputsInfo.Name,LinearModel.OutputsInfo.OutputPosition);
MagnSensIdx = signalIndexByName(MagnSensNames,LinearModel.OutputsInfo.Name,LinearModel.OutputsInfo.OutputPosition);
FluxCntrlIdx = signalIndexByName(FluxCntrlNames,LinearModel.OutputsInfo.Name,LinearModel.OutputsInfo.OutputPosition);


% Define gapIndex and strikeIndex, that is I'm removing all the NaN and
% negative gaps and NaN strikes
gapIndex       = [];
strikeIndex    = [];
FluxSensIndex  = [];
MagnSensIndex  = [];
FluxCntrlIndex = [];

XPntIndex      = [];

for i = 1:length(gapIdx)
    if ~isnan(LinearModel.YEquil(gapIdx(i))) && LinearModel.YEquil(gapIdx(i)) > 0
        gapIndex(end+1) = i;
    end
end

for i = 1:length(strikeIdx)
    if ~isnan(LinearModel.YEquil(strikeIdx(i))) 
        strikeIndex(end+1) = i;
    end
end

for i = 1:length(FluxSensIdx)
    if ~isnan(LinearModel.YEquil(FluxSensIdx(i))) 
       FluxSensIndex(end+1) = i;
    end
end

for i = 1:length(MagnSensIdx)
    if ~isnan(LinearModel.YEquil(MagnSensIdx(i))) 
      MagnSensIndex(end+1) = i;
    end
end

for i = 1:length(FluxCntrlIdx)
    if ~isnan(LinearModel.YEquil(FluxCntrlIdx(i))) 
      FluxCntrlIndex(end+1) = i;
    end
end



% Equilibrium values for the gaps
equilGaps      = LinearModel.YEquil(gapIdx);
% Equilibrium values for the strike 
equilStrikes   = LinearModel.YEquil(strikeIdx);
% Equilibrium values for the strike 
equilFluxSens  = LinearModel.YEquil(FluxSensIdx);
equilMagnSens  = LinearModel.YEquil(MagnSensIdx);
equilFluxCntrl = LinearModel.YEquil(FluxCntrlIdx);

targetGaps = targetConf.LinearModel.YEquil(gapIdx);
for i = 1 : nOfGaps
    Rt(i, 1) = Input_struct.r_sens_gap(i) + targetGaps(i)*cosd(Input_struct.theta_sens_gap_deg(i));
    Zt(i, 1) = Input_struct.z_sens_gap(i) + targetGaps(i)*sind(Input_struct.theta_sens_gap_deg(i));
end

%% Initialize X-point
xIdx = signalIndexByName({'CV-RX', 'CV-ZX', 'CV-FLUX'}, LinearModel.OutputsInfo.Name, 1:length(LinearModel.OutputsInfo.Name));
 equilRxZx= LinearModel.YEquil(xIdx);
equilRxZx_sqzd=targetConf.LinearModel.YEquil(xIdx);

%%% for simulink reference
Reference=[    startTime-neg_time           ,zeros(1,nOfFluxCntrl-1), equilRxZx(1:2)'
                startTime                   ,zeros(1,nOfFluxCntrl-1) ,equilRxZx(1:2)'
                startTime+transitionTime    ,zeros(1,nOfFluxCntrl-1) ,equilRxZx_sqzd(1:2)'];
            

            
% Plots the first wall and the equilbrium 
close all
figure('Position', [0, 0, 500, 750])
pdemesh(Input_struct.p, Input_struct.e, []), hold on, axis([0, 6, -4.5, 4.5]), axis equal
pdecont(Input_struct.p,Input_struct.t,x_np(1:length(Input_struct.p)),y_np(strcmp(y_type,'psb_c'))*[1 1]);
xlabel('R[m]', 'Interpreter', 'latex')
ylabel('Z[m]', 'Interpreter', 'latex')
hold on

% target shape
H = pdecont(targetEquil.Input_struct.p,targetEquil.Input_struct.t,targetEquil.x_np(1:length(targetEquil.Input_struct.p))...
    ,targetEquil.y_np(strcmp(targetEquil.y_type,'psb_c'))*[1 1]);
set(H, 'Color', 'magenta')
set(H,'LineWidth',2)
set(H,'DisplayName','Limiter Plasma 3.5 [MA]')
hold on


%% Build state space model
VDNames = {'CS1','CS2','CS3','CS4','EF1','EF2','EF3','EF4','EF5','EF6','VSU','VSL','Vpl'}; 
PFNames = {'CS1','CS2','CS3','CS4','EF1','EF2','EF3','EF4','EF5','EF6','VSU','VSL'};
OutputNames = {PFNames{1:end},'Ipl',gapNames{1:end},strikeNames{1:end},'CV-RX','CV-ZX','CV-FLUX','CV-ZC',FluxSensNames{1:end},MagnSensNames{1:end},FluxCntrlNames{1:end}};
outputIdx = signalIndexByName(OutputNames,LinearModel.OutputsInfo.Name,LinearModel.OutputsInfo.OutputPosition);

A = -inv(LinearModel.L)*LinearModel.R;
B = inv(LinearModel.L);
E = -inv(LinearModel.L)*LinearModel.LE;
% [VV,DD] = eig(A);
% ii = find(diag(DD)>0);
% DD(ii,ii) = -DD(ii,ii);
% Astab = VV*DD*inv(VV);
% 
% A_contr = Astab;
A_contr = A;
B_contr = [B(:,[1:12 end]) A_contr*E];
C_contr = LinearModel.C(outputIdx,:);
D_contr = [LinearModel.D(outputIdx,[1:12 end]) LinearModel.F(outputIdx,:)+C_contr*E];


%% Initial state
xi0 = zeros(1, size(A_contr, 1)); % substitute a VDE-related x0 if needed

%% Initialize plasma current
Ip0 = LinearModel.XEquil(LinearModel.PlasmaCurrentInfo.StatePosition);
vLoop = LinearModel.R(end,end)*Ip0;
IpIdx = signalIndexByName({'Ipl'}, OutputNames, 1:length(OutputNames));

%% Equilibrium inputs and outputs 
equilOutputs = LinearModel.YEquil(outputIdx); % from model outputs  
equilDist = LinearModel.YEquil(LinearModel.DisturbancesInfo.OutputPosition);
equilVolts = zeros(13,1);

%% Build reference waveforms for Ip and shape and scenario currents
% Ip_ref
initialIp = Ip0;
finalIp = targetConf.LinearModel.XEquil(LinearModel.PlasmaCurrentInfo.StatePosition);

Ip_ref = [startTime-neg_time initialIp;  startTime initialIp;startTime+transitionTime finalIp];
% Ip_ref = [startTime-neg_time initialIp;  startTime initialIp;startTime+transitionTime finalIp; 5 5e6];

% Gap_ref
ind1=14;
ind2=14+84;
tempOutputs = targetConf.LinearModel.YEquil(outputIdx);


all_initialGaps = equilOutputs([ind1:ind2]);
all_finalGaps = tempOutputs([ind1:ind2]);


if nOfFluxCntrl == 9
% %% Japanese selection 1 (Add one more gap at the begining)
japs_gaps=[13,13,22,36,52,66,76,85,84]';
japs_gaps=japs_gaps+ind1-1;
initialGaps = equilOutputs(japs_gaps);
finalGaps = tempOutputs(japs_gaps);
cntrl_gaps=[13,13,22,36,52,66,76,84,85]';
NoControlledGaps='Controlled Gaps = 8'
%%%%% Obtain R,Z coordinates of the equilibrium gaps
ii=length(Input_struct.r_sens_gap);

r_gap=Input_struct.r_sens_gap([13,13,22,36,52,66,76,84,85]);
 r_gap_sqzd=targetEquil.Input_struct.r_sens_gap([13,13,22,36,52,66,76,84,85]);
    
z_gap=Input_struct.z_sens_gap([13,13,22,36,52,66,76,84,85]);
z_gap_sqzd=targetEquil.Input_struct.z_sens_gap([13,13,22,36,52,66,76,84,85]);

theta_cntrl_gap=Input_struct.theta_sens_gap_deg([13,13,22,36,52,66,76,84,85]);
theta_cntrl_gap_sqzd=targetEquil.Input_struct.theta_sens_gap_deg([13,13,22,36,52,66,76,84,85]);

 elseif nOfFluxCntrl ==7


% % %Japanese selection 2
japs_gaps=[21,21,36,47,66,84,85]';
japs_gaps=japs_gaps+ind1-1;
initialGaps = equilOutputs(japs_gaps);
finalGaps = tempOutputs(japs_gaps);
cntrl_gaps=[21,21,36,47,66,84,85]';
NoControlledGaps='Controlled Gaps = 6'
%%%%% Obtain R,Z coordinates of the equilibrium gaps
ii=length(Input_struct.r_sens_gap);

r_gap=Input_struct.r_sens_gap([21,21,36,47,66,84,85]);
r_gap_sqzd=targetEquil.Input_struct.r_sens_gap([21,21,36,47,66,84,85]);

z_gap=Input_struct.z_sens_gap([21,21,36,47,66,84,85]);
z_gap_sqzd=targetEquil.Input_struct.z_sens_gap([21,21,36,47,66,84,85]);

theta_cntrl_gap=Input_struct.theta_sens_gap_deg([21,21,36,47,66,84,85]);
 theta_cntrl_gap_sqzd=targetEquil.Input_struct.theta_sens_gap_deg([21,21,36,47,66,84,85]);
else


%%Selection of 20 gaps
initialGaps = equilOutputs([round((ind1+1:4.6:(ind2-2))),ind2-1,ind2]);
finalGaps = tempOutputs([round((ind1+1:4.6:(ind2-2))),ind2-1,ind2]);
cntrl_gaps=round((ind1+1:4.6:(ind2-2)))-13;
NoControlledGaps='Controlled Gaps = 20'

%%%%% Obtain R,Z coordinates of the equilibrium gaps
ii=length(Input_struct.r_sens_gap);
r_gap=Input_struct.r_sens_gap([round((1+1:4.6:(ii-10))),ii-9,ii-8]);
r_gap_sqzd=targetEquil.Input_struct.r_sens_gap([round((1+1:4.6:(ii-10))),ii-9,ii-8]);

z_gap=Input_struct.z_sens_gap([round((1+1:4.6:(ii-10))),ii-9,ii-8]);
z_gap_sqzd=targetEquil.Input_struct.z_sens_gap([round((1+1:4.6:(ii-10))),ii-9,ii-8]);

theta_cntrl_gap=Input_struct.theta_sens_gap_deg([round((1+1:4.6:(ii-10))),ii-9,ii-8]);
theta_cntrl_gap_sqzd=targetEquil.Input_struct.theta_sens_gap_deg([round((1+1:4.6:(ii-10))),ii-9,ii-8]);

%%
r_cntrl_gap_eq=r_gap+initialGaps.*cosd(theta_cntrl_gap);
z_cntrl_gap_eq=z_gap+initialGaps.*sind(theta_cntrl_gap);

r_cntrl_gap_eq_sqzd=r_gap_sqzd+finalGaps.*cosd(theta_cntrl_gap_sqzd);
z_cntrl_gap_eq_sqzd=z_gap_sqzd+finalGaps.*sind(theta_cntrl_gap_sqzd);
end



%% Reference of the controil points for CCS

r_cntrl_gap_eq=r_gap+initialGaps.*cosd(theta_cntrl_gap);
z_cntrl_gap_eq=z_gap+initialGaps.*sind(theta_cntrl_gap);

r_cntrl_gap_eq_sqzd=r_gap_sqzd+finalGaps.*cosd(theta_cntrl_gap_sqzd);
z_cntrl_gap_eq_sqzd=z_gap_sqzd+finalGaps.*sind(theta_cntrl_gap_sqzd);



Gap_ref = [startTime-neg_time initialGaps'; startTime initialGaps';startTime+transitionTime finalGaps'];
r_cntrl_ref=[startTime-neg_time r_cntrl_gap_eq'; startTime r_cntrl_gap_eq';startTime+transitionTime r_cntrl_gap_eq_sqzd'];
z_cntrl_ref=[startTime-neg_time z_cntrl_gap_eq'; startTime z_cntrl_gap_eq';startTime+transitionTime z_cntrl_gap_eq_sqzd'];



%% I_scenario
initialPFCurr = equilOutputs([1:10]);
finalPFCurr = tempOutputs([1:10]);
IPF_scenario = [startTime-neg_time initialPFCurr';startTime initialPFCurr';startTime+transitionTime finalPFCurr'];


% Define the disturbance waveforms

%% Disruptions
load('../../Data from the Japanese/Temporal_evolutions_of_betap_and_li/li_bp_IpRU2.mat');
load('../../disturb_PID');
%betap = [0 0*Ip0+equilDist(1); 0.1 0*Ip0+equilDist(1); 0.1+1e-3 -0.14*Ip0+equilDist(1); 0.2 -0.14*Ip0+equilDist(1)];
%li = [0 0*Ip0+equilDist(2); 0.1 0*Ip0+equilDist(2); 0.1+3e-3 -0.15*Ip0+equilDist(2); 0.2 -0.15*Ip0+equilDist(2)];

betap_ur=betap_Ur2+equilDist(1)/Ip0;
li_ur=li_Ur2+equilDist(2)/Ip0;

betap_miy=betap_Miy2+equilDist(1)/Ip0;
li_miy=li_Miy2+equilDist(2)/Ip0;

betap_ur=Ip0*betap_ur;
li_ur=Ip0*li_ur;

 li_ur=[Ura_t2;li_ur]';
 betap_ur=[Ura_t2;betap_ur]';
 
 
betap_miy=Ip0*betap_miy;
li_miy=Ip0*li_miy;
li_miy=[Miy_t2;li_miy]';
betap_miy=[Miy_t2;betap_miy]';
% 
% betap=betap_miy;
% li=li_miy;
% type_dis='Miyata'

if dist==1
    
betap=betap_ur;
li=li_ur;
type_dis='Urano'

betap = [startTime-3 betap(1,2);betap];
li = [startTime-3 li(1,2);li];

elseif dist==2
    
 compELM.betap(:,2)=Ip0*compELM.betap(:,2);%
 compELM.li(:,2)=Ip0*compELM.li(:,2);
 betap = [startTime-3 equilDist(1);startTime equilDist(1);compELM.betap];
 li = [startTime-3 equilDist(2);startTime equilDist(2);compELM.li];
type_dis='comp_ELM'

elseif dist==3

 ELM.betap(:,2)=Ip0*ELM.betap(:,2);%
 betap = [startTime-3 equilDist(1);startTime equilDist(1);ELM.betap];
 li = [startTime-3 equilDist(2);startTime equilDist(2)];
type_dis='ELM'
    
elseif dist==4
    
     clear betap li
 minor1.betap(:,2)=Ip0*minor1.betap(:,2);
 minor1.li(:,2)=Ip0*minor1.li(:,2);
 betap = [startTime-3 equilDist(1);startTime equilDist(1);minor1.betap];
 li = [startTime-3 equilDist(2);startTime equilDist(2);minor1.li];
 type_dis='minor_disr'

else
    
    %%%% New minor disruption definition ITER 
     type_dis='minor_disr_ITERlike'
clear betap li
minorITER.betap=[0 0; 0.1 0; 0.1+1e-3 -0.53*equilDist(1); 0.2 -0.53*equilDist(1); stopTime -0.53*equilDist(1)];
minorITER.betap(:,2) = equilDist(1)+minorITER.betap(:,2);
minorITER.li = [0 0; 0.1 0; 0.1+1e-3 -0.125*equilDist(2); 0.2 -0.125*equilDist(2); stopTime -0.125*equilDist(2)];
minorITER.li(:,2) = equilDist(2)+minorITER.li(:,2);

 betap = [startTime-3 equilDist(1);startTime equilDist(1);minorITER.betap];
 li = [startTime-3 equilDist(2);startTime equilDist(2);minorITER.li];

end


%% Controllers 

% Load the PF current controller
load('../../KPFcurr.mat','KPFcurr');

%% Model of the PF power supplies
% Set limits for the simulation scheme
% VSsat = 1e3;
CSsat = [1 1.3 1.3 1] * 1e3;
EFsat = [1 .97 .97 .97 .97 1] * 1e3;

I_sat_upper = [2 2 2 2 1 1 2 2 1 1]*1e4;
I_sat_lower = [-2 -2 -2 -2 -2 -2 -2 -2 -2 -2]*1e4; 

%% PFC control
PSdelay = 1.5e-3; 
PStau = 3e-3;
PSmodel = tf([1],[PStau 1]) * eye(nPF);

%% Ip control
% PF combination that produces the transformer field
kIp = -[7.7175
        6.6999
        6.8024
        7.6938
        0.3730
        0.4150
        4.5207
        2.9209
        0.7198
        0.2332]*1e3;  
kIp = kIp/norm(kIp);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PID controller
IpControl.Kp = 0.001;
IpControl.Ki = 100 * IpControl.Kp;
IpControl.Kd = 0;
IpControl.Td = 1;

Ip_PID = buildPID(IpControl);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% XSC matrix
% load('../../JT-60SA simulation scheme/XSC_SOF.mat','XSC_matrix')
% load('../../XSC_SOF.mat','XSC_matrix')
% load('XSC_SOF_flux.mat','XSC_matrix')
XSC_isoflux_FBC
% XSC PID controller
XSC_PID.Kp =0.05* ones(1, nPF);  
XSC_PID.Ki = 200*XSC_PID.Kp;
XSC_PID.Kd = zeros(1, nPF);
XSC_PID.Td = ones(1, nPF); % dummy value

% % % XSC PID controller
% XSC_PID.Kp = 11 * ones(1, nPF);  
% XSC_PID.Ki = 450 * ones(1, nPF);
% XSC_PID.Kd = 29*ones(1, nPF);
% XSC_PID.Td = (1/2000)*ones(1, nPF); % dummy value


XSC_PI = buildPID(XSC_PID);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% VS system
VSL_Vsat = [-1, +1] * 1e3; % ************ In-vessel coils
VSU_Vsat = [-1, +1] * 1e3;  
VSL_Isat = [-5, +5] * 1e3; 
VSU_Isat = [-5, +5] * 1e3; 

VS.ICGain = -0.0076; % current gain
VS.ZdotGain = 20/5.5e6*Ip0; % zdot gain
VS.OverallGain = -40;
VS.Filter = tf([1, 0], [1/1000, 1]);


return





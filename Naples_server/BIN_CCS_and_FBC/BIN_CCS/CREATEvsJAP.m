system('./run.csh')
t=0;
Ip=equilOutputs(13);
IPF=equilOutputs(1:10);
IIC=equilOutputs(11:12);
psi_sens=-equilOutputs(110:143)/2/pi;
B_sens=equilOutputs(144:188);
Ip = Ip*1e-6;
IPF = IPF*1e-6;
IIC = IIC*1e-6;


%%% Control points for isoflux
ncp = 6;
% rcp = [0 1:10]; % 1st value must be zero?
rcp=[0.00000000e+00, 4.00000000e+00, 3.60000000e+00, 2.10000000e+00 ,...
2.10000000e+00, 3.60000000e+00];
% zcp = [0 11:20];
zcp=[ 0.00000000e+00, 0.00000000e+00, 1.10000000e+00,...
1.60000000e+00,-1.60000000e+00,-1.10000000e+00];


% run CCS with CREATE equilibrium fluxes and fields
matlab2namelist_smlnk2(t,Ip,IPF,IIC,B_sens,psi_sens,ncp,rcp,zcp);
status=system('./ccs_sa');

% % %%% Field probes
% psi_sens1=[-2.69113000e-01, 1.85095000e-01, 1.33093100e+00, 1.87941100e+00, ...
%             2.09069200e+00, 2.16329800e+00, 2.08715300e+00, 1.95969000e+00, ...
%             2.35978800e+00, 2.24130500e+00, 1.77507700e+00, 1.19047200e+00, ...
%            -7.17410000e-02,-2.84081000e-01, 1.91600000e-02, 2.14502000e-01, ...
%             3.26301000e-01, 3.43940000e-01, 3.09876000e-01, 3.12354000e-01, ...
%             3.39005000e-01, 3.99300000e-01, 3.82358000e-01, 2.97814000e-01, ...
%            -5.46990000e-02,-4.07333000e-01,-4.07221000e-01, 1.03628300e+00, ...
%             6.39244000e-01, 5.96866000e-01, 4.45804000e-01, 6.30205000e-01, ...
%             6.51902000e-01, 5.26055000e-01];
% % 
% 
% %%% Flux probes
% B_sens1=[5.78728000e-01, 6.71539000e-01, 6.64365000e-01, 6.34128000e-01, ...
%      6.50877000e-01, 6.72249000e-01, 5.14968000e-01, 1.82327000e-01, ...
%     -3.20417000e-01,-3.79086000e-01,-8.64340000e-02, 1.59264000e-01, ...
%      4.55573000e-01, 7.63706000e-01, 1.05066000e+00, 1.12595100e+00, ...
%      1.06465700e+00, 7.86789000e-01, 4.86599000e-01, 2.05705000e-01,...
%     -3.01980000e-02,-3.03915000e-01,-3.37847000e-01,-5.59000000e-04, ...
%      4.12192000e-01, 6.97262000e-01, 8.46789000e-01, 7.63704000e-01, ...
%      6.47178000e-01, 5.49362000e-01, 5.95425000e-01, 6.70360000e-01, ...
%      7.15418000e-01, 7.66295000e-01, 7.89732000e-01, 7.83530000e-01, ...
%      7.46713000e-01, 7.62091000e-01, 7.06427000e-01, 7.49722000e-01, ...
%      7.75857000e-01, 7.76302000e-01, 7.48523000e-01, 7.00087000e-01, 6.15273000e-01];


% % check flux
% figure
% bar(1:34, [psi_sens psi_sens1'])
% title('flux')

% % check field
% figure
% bar(1:45, [B_sens B_sens1'])
% title('field')


%% Plot results
close all
figSens = figure('Position', [0, 0, 500, 750]);

% plot flux sensors
figure(figSens)
iii = find(isnan(Input_struct.theta_sens));
rr = Input_struct.r_sens(iii);
zz = Input_struct.z_sens(iii);
hold on
plot(rr,zz,'sr')

% plot field sensors
figure(figSens)
iii = find(not(isnan(Input_struct.theta_sens)));
rr = Input_struct.r_sens(iii);
zz = Input_struct.z_sens(iii);
tt = Input_struct.theta_sens(iii);
% Lets put it other direction like the japanese want
tt = Input_struct.theta_sens(iii)+pi;
hold on
quiver(rr,zz,.05*cos(tt),.05*sin(tt))

% target shape
hold on
H1 = pdecont(targetEquil.Input_struct.p,targetEquil.Input_struct.t,targetEquil.x_np(1:length(targetEquil.Input_struct.p))...
    ,targetEquil.y_np(strcmp(targetEquil.y_type,'psb_c'))*[1 1]);

H1(2).Visible = 'off';
delete(H1(2));
H1 = H1(1);

set(H1, 'Color', 'b')
set(H1,'LineWidth',2)
hold on



% plot CCS boundary
hold on
run('./m1.m')

% legend
%legend('\psi probes', 'B probes', 'CREATE LCFS', 'CCS LCFS')

% plot mesh
hold on
h = pdemesh(Input_struct.p, Input_struct.e, []), hold on, axis([0, 6, -4.5, 4.5]), axis equal;
set(h,'Color','k')


xlabel('R[m]', 'Interpreter', 'latex')
ylabel('Z[m]', 'Interpreter', 'latex')
axis equal
xlim([1 6])
ylim([-4 4])









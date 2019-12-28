% Matlab to fortran interface: write magnetic data into Fortran format
%--------------------------------------------------------------------------
%% Sample namelist
% &namelist70
% ctimeinp= 2.50000000e+01,
% aipinp= 5.49922000e+00,
% ccurinp=-2.30280000e-03,-1.01999000e-02,-1.17252000e-02,-3.15530000e-03,-1.06071000e-02,-1.39031000e-02, 1.34648000e-02, 1.08751000e-02, 3.66150000e-03,-1.91788000e-02, 0.00000000e+00, 0.00000000e+00,
% tprbinp= 5.98952000e-01, 6.97699000e-01, 6.87706000e-01, 6.52043000e-01, 6.63339000e-01, 6.80547000e-01, 5.21935000e-01, 1.90353000e-01,-3.08526000e-01,-3.58692000e-01,-6.10370000e-02, 1.78166000e-01, 4.46470000e-01, 7.18590000e-01, 9.68875000e-01, 1.03974400e+00, 9.97255000e-01, 7.62042000e-01, 4.98535000e-01, 2.52813000e-01,-1.56610000e-02,-3.57898000e-01,-3.96900000e-01, 3.70470000e-02, 4.64171000e-01, 7.29992000e-01, 8.63751000e-01, 7.76212000e-01, 6.59644000e-01, 5.64770000e-01, 6.11165000e-01, 6.88289000e-01, 7.35463000e-01, 7.88178000e-01, 8.13143000e-01, 8.07897000e-01, 7.70869000e-01, 7.86242000e-01, 7.27307000e-01, 7.69434000e-01, 7.95327000e-01, 7.95337000e-01, 7.68256000e-01, 7.23220000e-01, 6.41738000e-01,
% fluxlinp=-3.04888000e-01, 1.43377000e-01, 1.27811100e+00, 1.82558600e+00, 2.04106500e+00, 2.12951900e+00, 2.07320600e+00, 1.94239900e+00, 2.32280800e+00, 2.23463300e+00, 1.78696500e+00, 1.20675500e+00,-1.72498000e-01,-3.16095000e-01,-2.31020000e-02, 1.37910000e-01, 2.08354000e-01, 2.07783000e-01, 1.78870000e-01, 1.92262000e-01, 2.23462000e-01, 2.89419000e-01, 2.84047000e-01, 2.14084000e-01,-1.58541000e-01,-5.85912000e-01,-5.84641000e-01, 9.71039000e-01, 5.60386000e-01, 5.17306000e-01, 3.67702000e-01, 5.65279000e-01, 5.82352000e-01, 4.57933000e-01,
% nfluxinp=6,
% rfluxinp= 0.00000000e+00, 4.00000000e+00, 3.60000000e+00, 2.10000000e+00, 2.10000000e+00, 3.60000000e+00,
% zfluxinp= 0.00000000e+00, 0.00000000e+00, 1.10000000e+00, 1.60000000e+00,-1.60000000e+00,-1.10000000e+00,
% /

%% Explanation (name, number of slots, explanation, measurement unit)
% ctimeinp (1)   = computation time [s]
% aipinp   (1)   = plasma current [MA]
% ccurinp  (20)  = measured coil currents; 1-10--> PFC, 11-12-->FPPC [MA]
% tprbinp  (300) = field probes [T]
% fluxlinp (300) = flux loops [Wb/rad]
% nfluxinp (1)   = number of control points including X or limiter points [-]
% rfluxinp (20)  = R coordinate of control points, including X or limiter [m] (the 1st one is the X or limiter point, in the sample file it is dummy=0)
% zfluxinp (20)  = Z coordinate of control points, including X or limiter [m] (the 1st one is the X or limiter point, in the sample file it is dummy=0)
%%% NOTA: rfluxinp e zfluxinp sono da verificare (dice che vengono da
%%% rfluxrfm e zfluxrfm di fort.21, che � un file di FBC)
%--------------------------------------------------------------------------

%% Dummy data
t = 0;
Ip = 5.5e6;
IPF = [1 2 3 4 5 6 7 8 9 10]*1e3;
IIC = [11 12]*1e3;
B_sens = linspace(1e-1,5e-1,45); % 45 field sensors
psi_sens = linspace(-0.5,2.5,34); % 34 flux loops
ncp = 10+1;
rcp = [0 1:10]; % 1st value must be zero?
zcp = [0 11:20];


%% Write file
filename = 'pippo.70';

% Open file
fileID = fopen(filename,'w');
fprintf(fileID,'&namelist70\n');

% Time instant
fprintf(fileID,'ctimeinp=%e,\n',t);

% Plasma Current
fprintf(fileID,'aipinp=%e,\n',Ip);

% PF Currents
fprintf(fileID,'ccurinp=');
for i = 1 : length(IPF)
    fprintf(fileID,'%e,',IPF(i));
end
for i = 1 : length(IIC)
    fprintf(fileID,'%e,',IIC(i));
end
fprintf(fileID,'\n')

% Field sensors
fprintf(fileID,'fluxlinp=');
for i = 1 : length(B_sens)
    fprintf(fileID,'%e,',B_sens(i));
end
fprintf(fileID,'\n')

% Flux Sensors
fprintf(fileID,'tprbinp=');
for i = 1 : length(psi_sens)
    fprintf(fileID,'%e,',psi_sens(i));
end
fprintf(fileID,'\n')

% Control points
fprintf(fileID,'nfluxinp=%e,\n',ncp);

fprintf(fileID,'rfluxinp=');
for i = 1 : length(rcp)
    fprintf(fileID,'%e,',rcp(i)); % First value is 0
end
fprintf(fileID,'\n')

fprintf(fileID,'zfluxinp=');
for i = 1 : length(zcp)
    fprintf(fileID,'%e,',zcp(i));
end
fprintf(fileID,'\n')

% Close file
fprintf(fileID, '/');
fclose(fileID);

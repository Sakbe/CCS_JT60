function matlab2namelist_smlnk(t,Ip,IPF,IIC,B_sens,psi_sens)

%%% Dummy data
% t = 0;
% Ip = 5.5e6;
% IPF = [1 2 3 4 5 6 7 8 9 10]*1e3;
% IIC = [11 12]*1e3;
% B_sens = linspace(1e-1,5e-1,45); % 45 field sensors
% psi_sens = linspace(-0.5,2.5,34); % 34 flux loops
ncp = 10+1;
rcp = [0 1:10]; % 1st value must be zero?
zcp = [0 11:20];


%% Write file
folder = '/home/domenica/Scrivania/japanese_stuff/BIN_CCS_and_FBC/BIN_CCS/FORT70/';
filename = [folder 'test.70'];

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

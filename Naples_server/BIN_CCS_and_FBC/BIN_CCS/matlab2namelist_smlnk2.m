function matlab2namelist_smlnk2(t,Ip,IPF,IIC,B_sens,psi_sens,ncp,rcp,zcp)

%% Write file
filename = 'fort.70';

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

% Field Sensors
fprintf(fileID,'tprbinp=');
for i = 1 : length(B_sens)
    fprintf(fileID,'%e,',B_sens(i));
end
fprintf(fileID,'\n')

% Flux sensors
fprintf(fileID,'fluxlinp=');
for i = 1 : length(psi_sens)
    fprintf(fileID,'%e,',psi_sens(i));
end
fprintf(fileID,'\n')



% Control points
fprintf(fileID,'nfluxinp=%d,\n',int8(ncp));

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
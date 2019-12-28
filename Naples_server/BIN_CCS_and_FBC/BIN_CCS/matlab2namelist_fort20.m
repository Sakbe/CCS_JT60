function matlab2namelist_fort20(t,Ip,IPF,IIC,ncp,cntrlTime,psb_jap,r_cntrl_pnts,z_cntrl_pnts,r_cent,z_cent,majorR,li)

%% Write file
filename = 'fort.20';

% Open file
fileID = fopen(filename,'w');
fprintf(fileID,'&namelist20\n');

% Time instant
fprintf(fileID,'t1=%e,\n',t);

%Control cycle
fprintf(fileID,'ctldelt=%e,\n',cntrlTime)

%Number of control point
fprintf(fileID,'nfluxrfm=%d,\n',int8(ncp))

%Poloidal magnetic flux at control point
fprintf(fileID,'phiconrfm=');
for i = 1 : length(psb_jap)
    fprintf(fileID,'%e,',psb_jap(i));
end
fprintf(fileID,'\n')

%R-coordinates of control points
fprintf(fileID,'rfluxrfm=');
for i = 1 : length(r_cntrl_pnts)
    fprintf(fileID,'%e,',r_cntrl_pnts(i));
end
fprintf(fileID,'\n')

% Z-coordinates of control points.
fprintf(fileID,'zfluxrfm=');
for i = 1 : length(z_cntrl_pnts)
    fprintf(fileID,'%e,',z_cntrl_pnts(i));
end
fprintf(fileID,'\n')

%Coil currents used for equilibrium calculation
fprintf(fileID,'ccurrfm=');
for i = 1 : length(IPF)
    fprintf(fileID,'%e,',IPF(i));
end
for i = 1 : length(IIC)
    fprintf(fileID,'%e,',IIC(i));
end
fprintf(fileID,'\n')

%R-coordinate of plasma current center
fprintf(fileID,'rjrfm= %e,\n',r_cent)

%Z-coordinate of plasma current center
fprintf(fileID,'zjrfm= %e,\n',z_cent)

% Plasma Current
fprintf(fileID,'aiprfm=%e,\n',Ip);

%Plasma major radius
fprintf(fileID,'rprfm=%e,\n',majorR);

%Internal inductace
fprintf(fileID,'sli=%e,\n',li);



% Close file
fprintf(fileID, '/');
fclose(fileID);
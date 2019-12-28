function [IPFc,VPFc]=getIVPFfort21()

IPFc = zeros(110,1);
VPFc=zeros(110,1);


coder.extrinsic('fscanf');
coder.extrinsic('regexp');
coder.extrinsic('eval');

fid = fopen('fort.21');
A = fscanf(fid,'%s');
fclose(fid);

i1=double(0);
i2=double(0);
i3=double(0);


i1 = regexp(A, 'ai=');
i2 = regexp(A, 'vcom=');
i3 = regexp(A, 'nfluxrfm=');

%check before 19 and 21
ai= A(i1+18:i2-1);
vcom=A(i2+20:i3-1);

   
 IPFc=removeComasFBC(ai);
 VPFc=removeComasFBC(vcom);


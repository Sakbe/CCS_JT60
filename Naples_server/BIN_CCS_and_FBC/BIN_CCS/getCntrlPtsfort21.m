function [rcp,zcp]=getCntrlPtsfort21()

rcp = zeros(20,1);
zcp=zeros(20,1);


coder.extrinsic('fscanf');
coder.extrinsic('regexp');
coder.extrinsic('eval');

fid = fopen('fort.21');
A = fscanf(fid,'%s');
fclose(fid);

i1=double(0);
i2=double(0);
i3=double(0);


i1 = regexp(A, 'rfluxrfm=');
i2 = regexp(A, 'zfluxrfm=');
i3 = regexp(A, 'aiprfmbfr=');

%check before 19 and 21
rcp_temp= A(i1+9:i2-1);
zcp_temp=A(i2+9:i3-1);

   
 rcp=removeComas(rcp_temp);
 zcp=removeComas(zcp_temp);


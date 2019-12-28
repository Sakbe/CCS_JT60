function[rjr,zjr]=getCentroidfort71()

rjr=0;
zjr=0;


coder.extrinsic('fscanf');
coder.extrinsic('regexp');
coder.extrinsic('eval');

fid = fopen('fort.71');
A = fscanf(fid,'%s');
fclose(fid);

i1=double(0);
i2=double(0);
i3=double(0);

i1 = regexp(A,'rjrfm=');
i2 = regexp(A,'zjrfm=');
i3=regexp(A,'eddyrfm=')

rjr_temp = A(i1+6:i2-2);
zjr_temp=A(i2+6:i3-2);

rjr =str2double(rjr_temp);
zjr =str2double(zjr_temp);


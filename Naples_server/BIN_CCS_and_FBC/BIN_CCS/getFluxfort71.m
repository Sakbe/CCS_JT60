
function [psib,psicntrl,r_cntrl_pnts,z_cntrl_pnts]=getFluxfort71()
psib = 0;
psicntrl = zeros(20,1);
r_cntrl_pnts=zeros(20,1);
z_cntrl_pnts=zeros(20,1);

coder.extrinsic('fscanf');
coder.extrinsic('regexp');
coder.extrinsic('eval');

fid = fopen('fort.71');
A = fscanf(fid,'%s');
fclose(fid);

i1=double(0);
i2=double(0);
i3=double(0);
i4=double(0);
i5=double(0);

indx=zeros(1,9);
% psicntrl_temp='';

i1 = regexp(A, 'phisurrfm=');
i2 = regexp(A, 'phiconrfm=');
i3 = regexp(A, 'aiprfm=');
i4 = regexp(A, 'rfluxrfm=');
i5 = regexp(A, 'zfluxrfm=');

psib_temp     = A(i1+10:i2-2);
psicntrl_temp = A(i2+10:i3-1); % the first should be the boundary flux
rfluxrfm_temp=A(i4+9:i5-1); %1st is the X-point
zfluxrfm_temp=A(i5+9:i1-1);


% eval(['psib_temp = ' psib_temp])
% eval(['psicntrl_temp = [' psicntrl_temp ']'])

 psib =str2double(psib_temp);
% psicntrl = double(psicntrl_temp(1:20))';

%    indx=double(regexp(psicntrl_temp,','));
   
   psicntrl=removeComas(psicntrl_temp);
   r_cntrl_pnts=removeComas(rfluxrfm_temp);
   z_cntrl_pnts=removeComas(zfluxrfm_temp);

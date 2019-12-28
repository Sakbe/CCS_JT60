function dataD=removeComas(dataS)
coder.extrinsic('regexp');
dataD=zeros(20,1);
indx=zeros(1,20);


indx=double(regexp(dataS,','));

ii = 1;
%format long;
for i = 1 : length(indx)
    dataD(i,1) = str2double(dataS(ii:indx(i)-1));
    ii = indx(i) + 1;
end

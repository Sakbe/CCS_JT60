function dataD=removeComas(dataS)
coder.extrinsic('regexp');
dataD=zeros(110,1);
indx=zeros(1,110);


indx=double(regexp(dataS,','));

ii = 1;
for i = 1 : length(indx)
    dataD(i,1) = str2double(dataS(ii:indx(i)-1));
    ii = indx(i) + 1;
end

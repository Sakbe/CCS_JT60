fid = fopen('fort.71');
A = fscanf(fid,'%s');
fclose(fid);

i1 = regexp(A, 'phisurrfm=');
i2 = regexp(A, 'phiconrfm=');
i3 = regexp(A, 'aiprfm=');

psib     = A(i1+11:i2-2);
psicntrl = A(i2+11:i3-2); % the first should be the boundary flux

eval(['psib = ' psib])
eval(['psicntrl = [' psicntrl ']'])
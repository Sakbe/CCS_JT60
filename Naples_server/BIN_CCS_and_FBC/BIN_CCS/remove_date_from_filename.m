dd = dir('./');
for i = 3 : length(dd)
    try
        str = dd(i).name(strfind(dd(i).name, ' ('):strfind(dd(i).name, 'UTC)')+3);
        movefile(dd(i).name, strrep(dd(i).name, str, ''));
    catch
    end
end
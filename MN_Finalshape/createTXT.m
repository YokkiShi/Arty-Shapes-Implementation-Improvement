myFile=dir('*.txt');
for i=1:length(myFile)
    name{i} = myFile(i).name(1:length(myFile(i).name));
    delete(name{i});
end
myFile=dir('*.circ');
for i=1:length(myFile)
    name{i} = myFile(i).name(1:length(myFile(i).name));
    delete(name{i});
end
myFile=dir('*.tri');
for i=1:length(myFile)
    name{i} = myFile(i).name(1:length(myFile(i).name));
    delete(name{i});
end
%myFile=dir('*.robustCH');
%for i=1:length(myFile)
    %name{i} = myFile(i).name(1:length(myFile(i).name));
    %delete(name{i});
%end
myFile=dir('*.rect');
for i=1:length(myFile)
    name{i} = myFile(i).name(1:length(myFile(i).name));
    delete(name{i});
end
myFile=dir('*.sup');
for i=1:length(myFile)
    name{i} = myFile(i).name(1:length(myFile(i).name));
    delete(name{i});
end
system('./sh2p');
%system('./sh_robustCHp');
system('./sh3'); 
for i = 1:length(allLeafNodesFiles) % for each file (1 pix file for 1 part)
    curFile = char(allLeafNodesFiles(i));
    circ=strcat(curFile(1:end-4), '.', 'circ');
    rect=strcat(curFile(1:end-4), '.', 'rect');
    tri=strcat(curFile(1:end-4), '.', 'tri');
    sup=strcat(curFile(1:end-4), '.', 'sup');
    robustCH=strcat(curFile(1:end-4), '.', 'robustCH');
    txt=strcat(curFile(1:end-4), '.', 'txt');
    string_circ='./dt_pixel_stats %s %s %s %s';
    str_circ=sprintf(string_circ,circ,curFile,'>>',txt);
    string_rect='./dt_pixel_stats %s %s %s %s';
    str_rect=sprintf(string_rect,rect,curFile,'>>',txt);
    string_tri='./dt_pixel_stats %s %s %s %s';
    str_tri=sprintf(string_tri,tri,curFile,'>>',txt);
    string_sup='./dt_pixel_stats %s %s %s %s';
    str_sup=sprintf(string_sup,sup,curFile,'>>',txt);
    string_robustCH='./dt_pixel_stats %s %s %s %s';
    str_robustCH=sprintf(string_robustCH,robustCH,curFile,'>>',txt);
    system(str_circ);
    system(str_rect);
    system(str_tri);
    system(str_sup);
    system(str_robustCH);

end
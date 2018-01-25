function a = readPIX(filename)
% function to read in a Nx2 list of vertices for shapes fitted by Paul's
% code.

a = [];
start=0;
tline='x';
fid = fopen(filename, 'r');

while 1
    
    tline = fgetl(fid);
    
    if isequal(filename(end-2:end),'sup')
        if ~isempty(strfind(tline, 'superellipse'))
            start = 1;
        end
        
    elseif isequal(filename(end-3:end),'circ')
        if ~isempty(strfind(tline, 'circle'))
            start = 1;
        end
            
    else
        if strfind(tline, 'list')
            start = 1;
            tline = fgetl(fid);
        end
    end
    
    % start recording the coordinate information
    if start == 1
        
        %special format for information of circle and superellipse
        if isequal(filename(end-3:end),'circ')
            a = sscanf(tline', '%*s %g %g %g %g', [1, inf]);
            break
        elseif isequal(filename(end-2:end),'sup')
            a = sscanf(tline', '%*s %g %g %g %g %g %g %g %g %g %g %g %g', [1, inf]);
            break
        else
            %[-1 -1] is the ending representation
            if ismember(str2num(tline),[-1 -1;-1 0],'rows')
            break;
        
        else
        a = [a; str2num(tline)];
        end
        end
        
    end
       
end

    
fclose(fid);
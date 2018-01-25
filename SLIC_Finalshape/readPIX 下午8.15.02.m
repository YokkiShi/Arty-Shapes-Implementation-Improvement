function a = readPIX(filename)
% function to read in a Nx2 list of vertices for shapes fitted by Paul's
% code.

a = [];
start=0;
tline='x';
fid = fopen(filename, 'r');

while 1
    
    tline = fgetl(fid);
    
        if strfind(tline, 'list')
            start = 1;
            tline = fgetl(fid);
        end
  
    
    % start recording the coordinate information
    if start == 1
        
        %special format for information of circle and superellipse
        
        
            %[-1 -1] is the ending representation
            if ismember(str2num(tline),[-1 -1;-1 0],'rows')
            break;
        
        else
        a = [a; str2num(tline)];
        end
        
        
    end
       
end

    
fclose(fid);
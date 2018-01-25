function writePIX(M, filename)
% function to write a binary mask into a .pix format file that Paul's shape
% fitting code operates on.

if nargin < 2
    filename = 'file.pix';
end

% get the boundary pixel locations of the mask
B = bwboundaries(M);
% pick the largest object!
val = 0;
for i = 1:length(B)
    if size(B{i},1) > val
        BB = B{i};
        val = size(B{i},1);
    end
end

fid = fopen(filename, 'wt');
fprintf(fid, '%s\n', 'pixel');
fprintf(fid, '%s\n', 'list: 1');
fprintf(fid, '   %u    %u\n', fliplr(BB)');
fprintf(fid, '%s\n', '-1 -1');
fclose(fid);

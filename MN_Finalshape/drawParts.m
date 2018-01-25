% function to display image/object partitions
function [MsLabelled, axesNum, rgbM_ALL, objectNodeNums] = drawParts(MsALL, objectNodeNums, segColours, frgb, axesNum)

Ms = MsALL(objectNodeNums);
% figure out overlaping areas, because two scale regions overlap
MsOverlap = zeros(size(Ms{1}));
for i = 1:length(Ms)
    MsOverlap = MsOverlap + Ms{i};
end
MsOverlap = MsOverlap > 1;
% figure; imshow(MsALL);
% DISCARD overlaping areas from each segmentation mask
for i = 1:length(Ms)
    Ms{i} = logical(Ms{i}-and(Ms{i},MsOverlap));
end
fr = frgb(:,:,1);
fg = frgb(:,:,2);
fb = frgb(:,:,3);
rgbM = [];
rgbMs = {};
% build a labelled matrix of segments
MsLabelled = zeros(size(Ms{1}));
for i = 1:length(Ms)
    
    curM = Ms{i};
    curMask = curM==1;
    % The boundary of the current segment. we want to the highlight the
    % boundary with a colour
    curBound = bwmorph(curMask, 'remove');
    
    rgbM(:,:,1) = logical(curMask-curBound) .* double(fr);
    rgbM(:,:,2) = logical(curMask-curBound) .* double(fg);
    rgbM(:,:,3) = logical(curMask-curBound) .* double(fb);
    curColour = segColours(i,:);
    
    rgbM(:,:,1) = rgbM(:,:,1) + curBound * curColour(1);
    rgbM(:,:,2) = rgbM(:,:,2) + curBound * curColour(2);
    rgbM(:,:,3) = rgbM(:,:,3) + curBound * curColour(3);
    %rgbM(:,:,1) = rgbM(:,:,1) +  curColour(:,:,1);
    %rgbM(:,:,2) = rgbM(:,:,2) +  curColour(:,:,2);
    %rgbM(:,:,3) = rgbM(:,:,3) +  curColour(:,:,3);
    %    h = figure(i); imshow( rgbM );axax
    rgbMs{end+1} = rgbM;
    MsLabelled = MsLabelled + objectNodeNums(i)*curM;
end
% Show all segments in the same image
rgbM_ALL = zeros(size(rgbM));
for i = 1:length(rgbMs)
    curM = rgbMs{i};
    rgbM_ALL(:,:,1) = rgbM_ALL(:,:,1) + curM(:,:,1);
    rgbM_ALL(:,:,2) = rgbM_ALL(:,:,2) + curM(:,:,2);
    rgbM_ALL(:,:,3) = rgbM_ALL(:,:,3) + curM(:,:,3);
end

axesNum=axes('Parent',figure);
if nargin > 4
    imshow(rgbM_ALL,'Parent',axesNum);
else
    figure; axesNum = axes; imshow(rgbM_ALL);
end

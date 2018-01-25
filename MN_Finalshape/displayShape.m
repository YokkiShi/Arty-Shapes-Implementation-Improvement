% Read in shape fitting results
%allLeafNodesShapes = {};
%pixelPosition={};

% Specify what shape to fit to each part
% From {'rect', 'hull', 'robustCH', 'skull3', 'tri', 'para', 'circ', 'sup'}
% Correponding parts: {reye, leye, rear, lear, mouth, face}
%partShapes = {'sup', 'sup', 'sup', 'sup', 'sup', ...
              %'sup', 'sup', 'sup', 'sup', 'sup', ...
              %'sup', 'sup', 'sup', 'sup', 'sup', ...
              %'sup', 'sup', 'sup', 'sup', 'sup', ...
              %'sup', 'sup', 'sup', 'sup', 'sup'};
%for i = 1:length(allLeafNodesFiles) % for each file (1 pix file for 1 part)
    %curFile = char(allLeafNodesFiles(i));
    %pixelPosition{end+1}=readPIX(strcat(curFile(1:end-4), '.', 'pix'));
    %allLeafNodesShapes{end+1} = readPIX(strcat(curFile(1:end-4), '.', partShapes{i}));
%end

% Overlay shapes fitted to the segmented object itself
drawParts(allLeafNodesMs, 1:length(allLeafNodesMs),segColours,frgb,2);
hold on;

imshow('i31.jpg');
for i = 1:length(allLeafNodesShapes)
    
    curShape = allLeafNodesShapes{i};
    switch partShapes{i}
        case {'rect', 'hull', 'robustCH', 'skull3', 'tri', 'para'}
            plot(curShape(:,1),curShape(:,2),'r-', 'LineWidth',3);
        case 'circ'
            drawCircle(curShape(2),curShape(3),curShape(4));
        case 'sup'
            drawSuperEllipse(curShape(2),curShape(3),curShape(8),curShape(9),curShape(10),curShape(11));           
    end
    hold on;
    
end
axis image;

% show only shapes fitted
figure;
for i = 1:length(allLeafNodesShapes)
    curShape = allLeafNodesShapes{i};
    switch partShapes{i}
        case {'rect', 'hull', 'robustCH', 'skull3', 'tri', 'para'}
            
            patch(curShape(:,1),curShape(:,2),segColours(i,:)/255);
            
        case 'circ'
            fillCircle(curShape(2),curShape(3),curShape(4),segColours(i,:)/255);
            
        case 'sup'
            fillSuperEllipse(curShape(2),curShape(3),curShape(8),curShape(9),curShape(10),curShape(11),segColours(i,:)/255);     
    end
    hold on;
end
axis image; axis off; axis ij;
%save('pixelPosition.mat','pixelPosition');
save('allLeafNodesFiles.mat','allLeafNodesFiles');
%save('partShapes.mat','partShapes');


partShapes={};
allLeafNodesShapes={};
for i = 1:length(allLeafNodesFiles)
    curFile = char(allLeafNodesFiles(i));
    curTxt=strcat(curFile(1:end-4), '.', 'txt');
    test=load(curTxt);
    if test(3,2)<= 0.031108
        class='tri';
    elseif test(1,4)>1.7115
        class='tri';
    elseif test(3,4)<=-1.7245
        class='tri';
    elseif test(2,2)>0.060578 && test(4,3)>-0.092613 && test(5,3)>1.8777
        class='tri';
    elseif test(5,3)<=0.3533
        class='tri';
    elseif test(2,4)>-1.2718 && test(3,3)<=0.16205 && test(4,3)>-0.092613 && test(5,3)>1.8777
        class='tri';
    elseif test(1,4) >0.02493 && test(5,1)<=0.011845
        class='tri';
    elseif test(2,2)>0.033323 && test(2,4)>-1.2718 && test(4,2)<=0.027773
        class='sup';
    elseif test(1,1)>0.18522 && test(4,3)>-0.092613 && test(4,3)<= 0.083263
        class='sup';
    elseif test(2,4)> -0.96067 && test(4,4)<= -1.8429
        class='sup';
    elseif test(1,2)<= 0.044915
        class='robustCH';
    elseif 0.05733<test(1,1)<=0.22306 &&test(2,2)>0.033323 &&test(4,1)<= 0.030966 &&test(4,2)> 0.027773
        class='rect';
    elseif 0.028309<test(2,2)<= 0.033323 && test(4,4)> -1.8429
         class='rect';
    elseif test(2,1)<= 0.03785 && test(2,2) > 0.033323
        class='rect';
    elseif 0.051541<test(1,1)<=0.05733 && test(4,4)> -1.8429
        class='rect';
    elseif test(2,4)<= -1.2718 && test(4,4)> -1.8429
        class='rect';
    elseif test(1,1)> 0.14752 && test(4,3)<= -0.092613
        class='rect';
    elseif test(1,4)<= -1.3336 && test(4,4)> -1.8429
        class='rect';
    elseif test(1,4)> 0.02493 && test(3,4)> -1.7245 && test(5,3) <= 0.69179
        class='rect';
    else
        class='robustCH';
    end
    
    partShapes{i}=class;
    allLeafNodesShapes{end+1}=readPIX(strcat(curFile(1:end-4), '.', class));
end


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

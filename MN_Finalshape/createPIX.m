myFile=dir('*.pix');
for i=1:length(myFile)
    name{i} = myFile(i).name(1:length(myFile(i).name));
    delete(name{i});
end
myFile=dir('*.tri');
for i=1:length(myFile)
    name{i} = myFile(i).name(1:length(myFile(i).name));
    delete(name{i});
end
myFile=dir('*.rect');
for i=1:length(myFile)
    name{i} = myFile(i).name(1:length(myFile(i).name));
    delete(name{i});
end
myFile=dir('*.circ');
for i=1:length(myFile)
    name{i} = myFile(i).name(1:length(myFile(i).name));
    delete(name{i});
end
myFile=dir('*.sup');
for i=1:length(myFile)
    name{i} = myFile(i).name(1:length(myFile(i).name));
    delete(name{i});
end
%myFile=dir('*.robustCH');
%for i=1:length(myFile)
    %name{i} = myFile(i).name(1:length(myFile(i).name));
    %delete(name{i});
%end
s=pwd;
cd ('../ncut-multiscale');
load ('after-seg.mat');
load ('nsegs.mat');
cd(s);

%create PIX files

allLeafNodesFiles={};
allLeafNodesMs={};%[matrix BW1,...,BW(nsegs)]
allLeafNodes={};

allPatchColours={};
segColours=[];
frgb=imread('i31.jpg');
fr = frgb(:,:,1);
fg = frgb(:,:,2);
fb = frgb(:,:,3);


for i=1:nsegs
    allLeafNodes{i}=(classes==i);
    eval(['BW',num2str(i),'=','(classes==i)',';']);
    curFileName = strcat('BW',num2str(i),'.pix');
    allLeafNodesFiles{end+1} = curFileName;
    writePIX(eval(['BW',num2str(i)]),curFileName);
    allLeafNodesMs{i}=eval(['BW',num2str(i)]);    
end 

for i=1:nsegs
    allPatchColours{i}(:,:,1)=allLeafNodesMs{i}.*double(fr);
    allPatchColours{i}(:,:,2)=allLeafNodesMs{i}.*double(fg);
    allPatchColours{i}(:,:,3)=allLeafNodesMs{i}.*double(fb);
end

for i=1:nsegs
    patchr=allPatchColours{i}(:,:,1);
    segColours(i,1)=mean(patchr(patchr~=0));
    patchg=allPatchColours{i}(:,:,2);
    segColours(i,2)=mean(patchg(patchg~=0));
    patchb=allPatchColours{i}(:,:,3);
    segColours(i,3)=mean(patchb(patchb~=0));
end

    %system('./sh1'); % parallelogram
    %system('./sh2'); % triangle, circle and rectangle
    %system('./sh3'); % superellipse
    %system('./sh_ch'); % convex hull
    %system('./sh_robustCH'); % robust convex hull (slow)
    %system('./sh_skull'); % convex skull (slow)


%structNodesNamesCell{face}
%allLeafNodesNames{reye,leye,rear..}
% Read in the shape fitting results!!




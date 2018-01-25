foreground=imread('foreground.jpg');
background=imread('background.jpg');
lab=rgb2lab(background);
l=lab(:,:,1); a=lab(:,:,2);b=lab(:,:,3);
lab1=rgb2lab(foreground);
l1=lab1(:,:,1); a1=lab1(:,:,2);b1=lab1(:,:,3);
E=sqrt((l-l1).^2+(a-a1).^2+(b-b1).^2)/4;
E(E<2.3)=0;
E1=logical(E);
foreground=double(foreground);
background=double(background);
foreground(:,:,1)=foreground(:,:,1).*E1;
foreground(:,:,2)=foreground(:,:,2).*E1;
foreground(:,:,3)=foreground(:,:,3).*E1;

test=imadd(foreground,background);
imshow(test/256);
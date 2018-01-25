function h=fillSuperEllipse(xcenter,ycenter,max,min,squareness,rotang,color)
hold on
t = -pi:pi/100:pi;
x = zeros(length(t),1);
y = x;

 x(:,1) = xcenter+sign(cos(t)).*cos(rotang).*max.*abs(cos(t)).^(2/squareness);
 y(:,1) = ycenter+sign(sin(t)).*sin(rotang).*min.*abs(sin(t)).^(2/squareness);

 h=fill(x,y,color);
 hold off
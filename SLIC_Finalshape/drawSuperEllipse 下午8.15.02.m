function h=drawSuperEllipse(xcenter,ycenter,max,min,squareness,rotang)
hold on
t = -pi:pi/100:pi;
x = zeros(length(t),1);
y = x;

 x(:,1) = xcenter+sign(cos(t)).*cos(rotang).*max.*abs(cos(t)).^(2/squareness);
 y(:,1) = ycenter+sign(sin(t)).*sin(rotang).*min.*abs(sin(t)).^(2/squareness);

 h=plot(x,y);
 hold off

%axis equal tight off

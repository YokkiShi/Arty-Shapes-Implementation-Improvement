function h=fillCircle(x,y,r,color)
hold on
th = 0:pi/50:2*pi;
xunit = r * cos(th) + x;
yunit = r * sin(th) + y;
h = fill(xunit, yunit,color);
hold off